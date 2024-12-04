# PFAPACK wrapper of the C library.

# This module wraps
# skpfa.o skpf10.o
# and not yet
# skbpfa.o skbpf10.o sktrf.o sktrd.o skbtrd.o

from __future__ import annotations

import ctypes
import os
from pathlib import Path
from typing import Final

import numpy as np
from numpy.ctypeslib import ndpointer

from pfapack.exceptions import (
    ComputationError,
    InvalidDimensionError,
    InvalidParameterError,
)


def _find_library() -> ctypes.CDLL:
    """Find and load the PFAPACK C library."""
    _folder: Final = Path(__file__).parent
    _build_folder: Final = _folder.parent / "build"

    # On Windows, load OpenBLAS first
    if os.name == "nt":
        openblas_path = _folder / "libopenblas.dll"
        if openblas_path.exists():
            try:
                ctypes.CDLL(str(openblas_path))
            except OSError as e:
                print(f"Warning: Failed to load OpenBLAS: {e}")
        else:
            print(f"Warning: OpenBLAS DLL not found at {openblas_path}")

    # Try all possible library names
    lib_names = [
        "cpfapack.dll",
        "libcpfapack.dll",
        "libcpfapack.so",
        "libcpfapack.dylib",
    ]

    # List of all possible paths
    possible_paths = []
    for lib_name in lib_names:
        possible_paths.append(_folder / lib_name)
        # Add build directories for editable install
        if _build_folder.exists():
            for p in _build_folder.glob("*"):
                if p.is_dir():
                    possible_paths.append(p / lib_name)

    # Try all possible paths
    errors = []
    for path in possible_paths:
        try:
            if path.exists():
                return ctypes.CDLL(str(path))
            else:
                errors.append(f"{path}: File does not exist")
        except OSError as e:
            errors.append(f"{path}: {e}")
            continue

    # If we get here, all attempts failed
    error_msg = "\n".join(
        [
            "Could not load PFAPACK library.",
            "Attempted paths:",
            *[f" {e}" for e in errors],
            f"Current directory: {os.getcwd()}",
            f"Package directory: {_folder}",
            f"Files in package directory: {list(_folder.glob('*'))}",
        ]
    )
    raise OSError(error_msg)


lib = _find_library()


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
) -> float | complex:
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

    Returns
    -------
    float | complex
        The Pfaffian of the matrix.

    Raises
    ------
    InvalidDimensionError
        If the matrix is not square or has odd dimensions.
    InvalidParameterError
        If uplo or method parameters are invalid.
    ComputationError
        If the computation fails.
    """
    if uplo not in ("U", "L"):
        raise InvalidParameterError(f"'uplo' must be 'U' or 'L', got {uplo!r}")
    if method not in ("P", "H"):
        raise InvalidParameterError(f"'method' must be 'P' or 'H', got {method!r}")

    uplo_bytes = uplo.encode()
    method_bytes = method.encode()

    # Check matrix is square
    if np.ndim(matrix) != 2 or np.shape(matrix)[0] != np.shape(matrix)[1]:
        raise InvalidDimensionError(
            f"Matrix must be square, got shape {np.shape(matrix)}"
        )

    # Check matrix has even dimensions
    n = np.shape(matrix)[0]
    if n % 2 != 0:
        raise InvalidDimensionError(f"Matrix dimension must be even, got {n}")

    if np.iscomplex(matrix).any():
        a = np.zeros((2,) + np.shape(matrix), dtype=np.float64, order="F")
        a[0] = np.real(matrix)
        a[1] = np.imag(matrix)
        if avoid_overflow:
            result_array = (ctypes.c_double * 4)(0.0, 0.0)
            success = skpf10_z(
                matrix.shape[0], a, result_array, uplo_bytes, method_bytes
            )
            x = result_array[0] + 1j * result_array[1]
            exp = result_array[2] + 1j * result_array[3]
            result = from_exp(x, exp)
        else:
            result_array = (ctypes.c_double * 2)(0.0, 0.0)
            success = skpfa_z(
                matrix.shape[0], a, result_array, uplo_bytes, method_bytes
            )
            result = result_array[0] + 1j * result_array[1]
    else:
        matrix_f = np.asarray(matrix, dtype=np.float64, order="F")
        if avoid_overflow:
            result_array = (ctypes.c_double * 2)(0.0, 0.0)
            success = skpf10_d(
                matrix.shape[0], matrix_f, result_array, uplo_bytes, method_bytes
            )
            result = from_exp(result_array[0], result_array[1])
        else:
            result_double = ctypes.c_double(0.0)
            success = skpfa_d(
                matrix.shape[0],
                matrix_f,
                ctypes.byref(result_double),
                uplo_bytes,
                method_bytes,
            )
            result = result_double.value

    if success != 0:
        raise ComputationError(f"PFAPACK returned error code {success}")
    return result
