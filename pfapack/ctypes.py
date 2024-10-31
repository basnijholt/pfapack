# PFAPACK wrapper of the C library.

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import importlib.resources as pkg_resources

try:
    # Accessing the resource in a way compatible with newer Python versions
    lib_resource = pkg_resources.files('pfapack').joinpath('libcpfapack.so')
    lib_path = lib_resource.resolve(strict=True)
    lib = ctypes.CDLL(str(lib_path))
except Exception as e:
    print(f"Error locating libcpfapack.so: {e}")
    raise

def _init(which):
    func = getattr(lib, which)
    func.restype = ctypes.c_int
    func.argtypes = [
        ctypes.c_int,
        ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    return func

def _init_batched(which):
    func = getattr(lib, which)
    func.restype = ctypes.c_int
    func.argtypes = [
        ctypes.c_int,  # batch_size
        ctypes.c_int,  # N
        ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),  # A_batch
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # PFAFF_batch
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    return func

def _init_batched_z(which):
    func = getattr(lib, which)
    func.restype = ctypes.c_int
    func.argtypes = [
        ctypes.c_int,  # batch_size
        ctypes.c_int,  # N
        ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),  # A_batch_real_imag
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # PFAFF_batch_real_imag
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    return func

def _init_batched_4d(which):
    func = getattr(lib, which)
    func.restype = ctypes.c_int
    func.argtypes = [
        ctypes.c_int,  # outer_batch_size
        ctypes.c_int,  # inner_batch_size
        ctypes.c_int,  # N
        ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),  # A_batch
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # PFAFF_batch
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    return func

def _init_batched_4d_z(which):
    func = getattr(lib, which)
    func.restype = ctypes.c_int
    func.argtypes = [
        ctypes.c_int,  # outer_batch_size
        ctypes.c_int,  # inner_batch_size
        ctypes.c_int,  # N
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # A_batch_real_imag - CHANGED THIS
        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),  # PFAFF_batch_real_imag
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    return func

functions = {
    "skpfa_d": _init("skpfa_d"),
    "skpf10_d": _init("skpf10_d"),
    "skpfa_z": _init("skpfa_z"),
    "skpf10_z": _init("skpf10_z"),
    "skpfa_batched_d": _init_batched("skpfa_batched_d"),
    "skpfa_batched_z": _init_batched_z("skpfa_batched_z"),
    "skpfa_batched_4d_d": _init_batched_4d("skpfa_batched_4d_d"),
    "skpfa_batched_4d_z": _init_batched_4d_z("skpfa_batched_4d_z"),
}

def from_exp(x, exp):
    try:
        return x * 10 ** exp
    except OverflowError:
        return x * np.inf

def pfaffian(matrices, uplo="U", method="P", avoid_overflow=False):
    uplo = uplo.encode()
    method = method.encode()
    if matrices.ndim < 2:
        raise ValueError("Input must be at least 2D.")

    N = matrices.shape[-1]
    if matrices.shape[-2] != N:
        raise ValueError("Last two dimensions must be square.")

    is_complex = np.iscomplexobj(matrices)
    dtype = np.float64 if not is_complex else np.complex128
    result = np.empty(matrices.shape[:-2], dtype=dtype)
    matrices = matrices.reshape(-1, N, N)
    func_suffix = '_z' if is_complex else '_d'

    for idx, matrix in enumerate(matrices):
        if is_complex:
            a = np.asfortranarray(matrix)
            a_complex = np.zeros((2,) + matrix.shape, dtype=np.float64, order="F")
            a_complex[0] = np.real(a)
            a_complex[1] = np.imag(a)
            a = a_complex
        else:
            a = np.asfortranarray(matrix, dtype=np.float64)

        if avoid_overflow:
            pfaffian_values = (ctypes.c_double * (4 if is_complex else 2))()
            success = functions[f"skpf10{func_suffix}"](N, a, pfaffian_values, uplo, method)
            if is_complex:
                x = pfaffian_values[0] + 1j * pfaffian_values[1]
                exp = pfaffian_values[2] + 1j * pfaffian_values[3]
                pfaffian = from_exp(x, exp)
            else:
                pfaffian = from_exp(pfaffian_values[0], pfaffian_values[1])
        else:
            pfaffian_value = (ctypes.c_double * (2 if is_complex else 1))()
            success = functions[f"skpfa{func_suffix}"](N, a, pfaffian_value, uplo, method)
            pfaffian = pfaffian_value[0] + (1j * pfaffian_value[1] if is_complex else 0)

        assert success == 0
        np.ravel(result)[idx] = pfaffian

    return result

def pfaffian_batched(matrices, uplo="U", method="P"):
    uplo = uplo.encode()
    method = method.encode()
    if matrices.ndim != 3:
        raise ValueError("Input must be 3D for batched operation.")

    batch_size, N, _ = matrices.shape
    if matrices.shape[-1] != N:
        raise ValueError("Last two dimensions of each matrix must be square.")

    is_complex = np.iscomplexobj(matrices)
    dtype = np.float64 if not is_complex else np.complex128
    result = np.empty(batch_size, dtype=dtype)

    if is_complex:
        # Ensure Fortran contiguous memory layout for complex data
        matrices = np.asfortranarray(matrices, dtype=np.complex128)
        matrices_c = np.empty((batch_size, 2, N, N), dtype=np.float64, order="F")
        matrices_c[:, 0, :, :] = np.real(matrices)
        matrices_c[:, 1, :, :] = np.imag(matrices)
        result_c = np.empty((batch_size, 2), dtype=np.float64)

        success = functions["skpfa_batched_z"](
            batch_size,
            N,
            matrices_c.ravel(),
            result_c.ravel(),
            uplo,
            method
        )
        result = result_c[:, 0] + 1j * result_c[:, 1]
    else:
        # Ensure Fortran contiguous memory layout for real data
        matrices = np.asfortranarray(matrices, dtype=np.float64)

        success = functions["skpfa_batched_d"](
            batch_size,
            N,
            matrices.ravel(),
            result,
            uplo,
            method
        )

    assert success == 0
    return result

def pfaffian_batched_4d(matrices, uplo="U", method="P"):
    uplo = uplo.encode()
    method = method.encode()
    if matrices.ndim != 4:
        raise ValueError("Input must be 4D for batched operation.")

    outer_batch_size, inner_batch_size, N, _ = matrices.shape
    if matrices.shape[-1] != N:
        raise ValueError("Last two dimensions of each matrix must be square.")

    is_complex = np.iscomplexobj(matrices)
    dtype = np.float64 if not is_complex else np.complex128
    result = np.empty((outer_batch_size, inner_batch_size), dtype=dtype)

    if is_complex:
        # Ensure Fortran contiguous memory layout for complex data
        matrices = np.asfortranarray(matrices, dtype=np.complex128)
        matrices_c = np.empty((outer_batch_size, inner_batch_size, 2, N, N), dtype=np.float64, order="F")
        matrices_c[:, :, 0, :, :] = np.real(matrices)
        matrices_c[:, :, 1, :, :] = np.imag(matrices)
        result_c = np.empty((outer_batch_size, inner_batch_size, 2), dtype=np.float64)

        success = functions["skpfa_batched_4d_z"](
            outer_batch_size,
            inner_batch_size,
            N,
            matrices_c.ravel(),
            result_c.ravel(),
            uplo,
            method
        )
        result = result_c[:, :, 0] + 1j * result_c[:, :, 1]
    else:
        # Ensure Fortran contiguous memory layout for real data
        matrices = np.asfortranarray(matrices, dtype=np.float64)

        success = functions["skpfa_batched_4d_d"](
            outer_batch_size,
            inner_batch_size,
            N,
            matrices.ravel(),
            result.ravel(),
            uplo,
            method
        )

    if success != 0:
        raise RuntimeError(f"PFAPACK returned error code {success}")

    return result

def pfaffian_batched_4d_cx(matrices, uplo="U", method="P"):
    uplo = uplo.encode()
    method = method.encode()
    if matrices.ndim != 5:
        raise ValueError("Input must be 5D for batched operation.")
    outer_batch_size, inner_batch_size, ncx, N, _ = matrices.shape
    if ncx != 2:
        raise ValueError("Unexpected layout of the input matrix.")
    if matrices.shape[-1] != N:
        raise ValueError("Last two dimensions of each matrix must be square.")

    result = np.empty((outer_batch_size, inner_batch_size), dtype=np.complex128)
    result_c = np.empty((outer_batch_size, inner_batch_size, 2), dtype=np.float64)

    success = functions["skpfa_batched_4d_z"](
        outer_batch_size,
        inner_batch_size,
        N,
        matrices,
        result_c,
        uplo,
        method
    )
    result = result_c[:, :, 0] + 1j * result_c[:, :, 1]
    if success != 0:
        raise RuntimeError(f"PFAPACK returned error code {success}")
    return result
