# PFAPACK wrapper of the C library.

# This module wraps
# skpfa.o skpf10.o
# and not yet
# skbpfa.o skbpf10.o sktrf.o sktrd.o skbtrd.o

import ctypes

import numpy as np
from numpy.ctypeslib import ndpointer

lib = ctypes.CDLL("libcpfapack.so")


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
        return x * 10**exp
    except OverflowError:
        return x * np.inf


def pfaffian(
    matrix: np.ndarray,
    uplo: str = "U",
    method: str = "P",
    avoid_overflow: bool = False,
):
    """Compute Pfaffian.

    Parameters
    ----------
    matrix : numpy.ndarray
        Square skew-symmetric matrix.
    uplo : str
        If 'U' ('L'), the upper (lower) triangle of the matrix is used.
    method : str
        If 'P' ('H'), the Parley-Reid (Householder) algorithm is used.
    avoid_overflow : bool
        If True, take special care to avoid numerical under- or
        overflow (at the cost of possible additional round-off errors).
    """
    uplo: bytes = uplo.encode()  # type: ignore[no-redef]
    method: bytes = method.encode()  # type: ignore[no-redef]
    assert np.ndim(matrix) == 2 and np.shape(matrix)[0] == np.shape(matrix)[1]
    if np.iscomplex(matrix).any():
        a = np.zeros((2,) + np.shape(matrix), dtype=np.float64, order="F")
        a[0] = np.real(matrix)
        a[1] = np.imag(matrix)
        if avoid_overflow:
            pfaffian = (ctypes.c_double * 4)(0.0, 0.0)
            success = skpf10_z(matrix.shape[0], a, pfaffian, uplo, method)
            x = pfaffian[0] + 1j * pfaffian[1]
            exp = pfaffian[2] + 1j * pfaffian[3]
            pfaffian = from_exp(x, exp)
        else:
            pfaffian = (ctypes.c_double * 2)(0.0, 0.0)
            success = skpfa_z(matrix.shape[0], a, pfaffian, uplo, method)
            pfaffian = pfaffian[0] + 1j * pfaffian[1]  # type: ignore[assignment]
    else:
        matrix = np.asarray(matrix, dtype=np.float64, order="F")
        if avoid_overflow:
            pfaffian = (ctypes.c_double * 2)(0.0, 0.0)
            success = skpf10_d(matrix.shape[0], matrix, pfaffian, uplo, method)
            pfaffian = from_exp(pfaffian[0], pfaffian[1])  # type: ignore[assignment]
        else:
            pfaffian = ctypes.c_double(0.0)  # type: ignore[assignment]
            success = skpfa_d(
                matrix.shape[0], matrix, ctypes.byref(pfaffian), uplo, method
            )
            pfaffian = pfaffian.value  # type: ignore[misc]
    assert success == 0
    return pfaffian
