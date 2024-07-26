# PFAPACK wrapper of the C library.

# This module wraps
# skpfa.o skpf10.o
# and not yet
# skbpfa.o skbpf10.o sktrf.o sktrd.o skbtrd.o

import ctypes

import numpy as np
from numpy.ctypeslib import ndpointer
import pkg_resources

# Try to find the library path using pkg_resources
try:
    lib_path = pkg_resources.resource_filename('pfapack', 'libcpfapack.so')
except Exception as e:
    print(f"Error locating libcpfapack.so: {e}")
    raise

# Load the library
lib = ctypes.CDLL(lib_path)

def _init(which):
    func = getattr(lib, which)
    func.restype = ctypes.c_int  # result type
    func.argtypes = [
        ctypes.c_int,
        ndpointer(ctypes.c_double, flags="F_CONTIGUOUS"),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    return func


skpfa_d = _init("skpfa_d")  # Pfaffian for real double
skpf10_d = _init("skpf10_d")
skpfa_z = _init("skpfa_z")  # Pfaffian for complex double
skpf10_z = _init("skpf10_z")


def from_exp(x, exp):
    """Convert pfapack overflow-safe representation (x, exponent) scalar number.

    Overflows are converted to infinities.
    """
    assert np.isclose(np.imag(exp), 0.0)
    try:
        return x * 10 ** exp
    except OverflowError:
        return x * np.inf


def pfaffian(matrices, uplo="U", method="P", avoid_overflow=False):
    """
    Compute Pfaffians for a batch of skew-symmetric matrices.

    Parameters
    ----------
    matrices : numpy.ndarray
        Batch of square skew-symmetric matrices, where computation is
        performed on the last two dimensions.
    uplo : str
        'U' for upper triangle, 'L' for lower triangle.
    method : str
        'P' for Parlett-Reid algorithm, 'H' for Householder.
    avoid_overflow : bool
        True to handle numerical under- or overflow.
    """
    uplo = uplo.encode()
    method = method.encode()

    if matrices.ndim < 2:
        raise ValueError("Input must be at least 2D.")

    N = matrices.shape[-1]
    if matrices.shape[-2] != N:
        raise ValueError("Last two dimensions must be square.")

    dtype = np.float64 if not np.iscomplexobj(matrices) else np.complex128
    result = np.empty(matrices.shape[:-2], dtype=dtype)
    matrices = matrices.reshape(-1, N, N)

    for idx, matrix in enumerate(matrices):
        if np.iscomplexobj(matrix):
            a = np.zeros((2,) + np.shape(matrix), dtype=np.float64, order="F")
            a[0] = np.real(matrix)
            a[1] = np.imag(matrix)
            if avoid_overflow:
                pfaffian_values = (ctypes.c_double * 4)(0.0, 0.0, 0.0, 0.0)
                success = skpf10_z(matrix.shape[0], a, pfaffian_values, uplo, method)
                x = pfaffian_values[0] + 1j * pfaffian_values[1]
                exp = pfaffian_values[2] + 1j * pfaffian_values[3]
                pfaffian = from_exp(x, exp)
            else:
                pfaffian_values = (ctypes.c_double * 2)(0.0, 0.0)
                success = skpfa_z(matrix.shape[0], a, pfaffian_values, uplo, method)
                pfaffian = pfaffian_values[0] + 1j * pfaffian_values[1]
        else:
            matrix = np.asarray(matrix, dtype=np.float64, order="F")
            if avoid_overflow:
                pfaffian_values = (ctypes.c_double * 2)(0.0, 0.0)
                success = skpf10_d(matrix.shape[0], matrix, pfaffian_values, uplo, method)
                pfaffian = from_exp(pfaffian_values[0], pfaffian_values[1])
            else:
                pfaffian_value = ctypes.c_double(0.0)
                success = skpfa_d(
                    matrix.shape[0], matrix, ctypes.byref(pfaffian_value), uplo, method
                )
                pfaffian = pfaffian_value.value
        
        assert success == 0
        np.ravel(result)[idx] = pfaffian

    return result

