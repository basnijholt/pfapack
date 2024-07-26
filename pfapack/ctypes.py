# PFAPACK wrapper of the C library.

# This module wraps
# skpfa.o skpf10.o
# and not yet
# skbpfa.o skbpf10.o sktrf.o sktrd.o skbtrd.o
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import pkg_resources

try:
    lib_path = pkg_resources.resource_filename('pfapack', 'libcpfapack.so')
    lib = ctypes.CDLL(lib_path)
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

functions = {
    "skpfa_d": _init("skpfa_d"),
    "skpf10_d": _init("skpf10_d"),
    "skpfa_z": _init("skpfa_z"),
    "skpf10_z": _init("skpf10_z")
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
