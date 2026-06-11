# PFAPACK wrapper of the C library.

# This module wraps
# skpfa.o skpf10.o
# and not yet
# skbpfa.o skbpf10.o sktrf.o sktrd.o skbtrd.o

from __future__ import annotations

import ctypes
import os
from concurrent.futures import ThreadPoolExecutor
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

    # On Windows, ensure OpenBLAS and its dependencies can be found
    if os.name == "nt":
        # Wheels repaired with delvewheel bundle OpenBLAS and the MinGW
        # runtime in pfapack.libs; no system installation is needed.
        _bundled = _folder.parent / "pfapack.libs"
        if _bundled.exists():
            os.add_dll_directory(str(_bundled))  # type: ignore[attr-defined]
            return _load_pfapack_dll(_folder, _build_folder)

        # Source builds: hunt for a system OpenBLAS.
        possible_blas_paths = [
            Path("C:/msys64/mingw64/bin"),  # MSYS2 MinGW64
            Path("C:/msys64/ucrt64/bin"),  # MSYS2 UCRT64
            Path("C:/msys64/clang64/bin"),  # MSYS2 Clang64
            Path(os.environ.get("OPENBLAS_PATH", "")).parent
            / "bin",  # Custom installation
            Path(os.environ.get("CONDA_PREFIX", "")) / "Library" / "bin",  # Conda
        ]

        # Try to find and load OpenBLAS from any of these locations
        blas_loaded = False
        for path in possible_blas_paths:
            if path.exists():
                try:
                    os.add_dll_directory(str(path))  # type: ignore[attr-defined]
                    openblas_path = path / "libopenblas.dll"
                    if openblas_path.exists():
                        ctypes.CDLL(str(openblas_path))
                        blas_loaded = True
                        break
                except OSError as e:
                    print(f"Warning: Failed to load OpenBLAS from {path}: {e}")

        if not blas_loaded:
            print("Warning: Could not load OpenBLAS from any known location")

    return _load_pfapack_dll(_folder, _build_folder)


def _load_pfapack_dll(_folder: Path, _build_folder: Path) -> ctypes.CDLL:
    """Load the PFAPACK C library from the package or build directory."""
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


def _init_batched(which, dtype):
    func = getattr(lib, which)
    func.restype = ctypes.c_int  # result type
    func.argtypes = [
        ctypes.c_int,  # batch_size
        ctypes.c_int,  # N
        ndpointer(dtype, flags="C_CONTIGUOUS"),  # A_batch
        ndpointer(dtype, flags="C_CONTIGUOUS"),  # PFAFF_batch
        ctypes.c_char_p,  # UPLO
        ctypes.c_char_p,  # MTHD
    ]
    return func


skpfa_d = _init("skpfa_d")  # Pfaffian for real double
skpf10_d = _init("skpf10_d")
skpfa_z = _init("skpfa_z")  # Pfaffian for complex double
skpf10_z = _init("skpf10_z")
skpfa_batched_d = _init_batched("skpfa_batched_d", np.float64)
skpfa_batched_z = _init_batched("skpfa_batched_z", np.complex128)


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


def pfaffian_batched(
    matrices: np.ndarray,
    uplo: str = "U",
    method: str = "P",
    workers: int = 1,
) -> np.ndarray:
    """Compute the Pfaffian of a batch of skew-symmetric matrices.

    Equivalent to calling :func:`pfaffian` on every matrix, but much
    faster for batches of many small matrices because the LAPACK
    workspace is allocated once for the whole batch.

    Ported from the batched implementation by Joshua Goings (@jjgoings),
    https://github.com/jjgoings/pfapack.

    Parameters
    ----------
    matrices : numpy.ndarray
        Array of shape ``(..., N, N)`` holding skew-symmetric matrices.
        Real input is computed in double precision, complex input in
        complex double precision; other dtypes are upcast.
    uplo : str
        If 'U' ('L'), the upper (lower) triangle of each matrix is used.
    method : str
        If 'P' ('H'), the Parley-Reid (Householder) algorithm is used.
    workers : int
        Number of threads to spread the batch over. The matrices are
        independent, so the speedup is roughly linear in the number of
        physical cores. -1 uses all available cores. Default is 1
        (serial); leave it at 1 if your application already
        parallelizes at a higher level.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(...,)`` with the Pfaffian of each matrix
        (a scalar if a single 2D matrix was passed). For matrices of
        odd dimension the Pfaffian is 0.

    Raises
    ------
    InvalidDimensionError
        If the input is not an array of square matrices.
    InvalidParameterError
        If uplo, method, or workers parameters are invalid.
    ComputationError
        If the computation fails.
    """
    if uplo not in ("U", "L"):
        raise InvalidParameterError(f"'uplo' must be 'U' or 'L', got {uplo!r}")
    if method not in ("P", "H"):
        raise InvalidParameterError(f"'method' must be 'P' or 'H', got {method!r}")
    if workers == -1:
        workers = os.cpu_count() or 1
    if not isinstance(workers, int) or workers < 1:
        raise InvalidParameterError(
            f"'workers' must be a positive integer or -1, got {workers!r}"
        )

    matrices = np.asarray(matrices)
    if matrices.ndim < 2 or matrices.shape[-1] != matrices.shape[-2]:
        raise InvalidDimensionError(
            "Expected an array of square matrices with shape (..., N, N), "
            f"got shape {matrices.shape}"
        )

    n = matrices.shape[-1]
    batch_shape = matrices.shape[:-2]
    batch_size = int(np.prod(batch_shape, dtype=np.int64))

    if np.iscomplexobj(matrices):
        dtype = np.complex128
        func = skpfa_batched_z
    else:
        dtype = np.float64
        func = skpfa_batched_d

    # Copies only if the dtype or memory layout requires it.
    a = np.ascontiguousarray(matrices, dtype=dtype).reshape(batch_size, n, n)
    pfaffians = np.empty(batch_size, dtype=dtype)
    uplo_bytes = uplo.encode()
    method_bytes = method.encode()

    def run(start: int, stop: int) -> int:
        return func(
            stop - start,
            n,
            a[start:stop],
            pfaffians[start:stop],
            uplo_bytes,
            method_bytes,
        )

    workers = min(workers, batch_size)
    if workers <= 1:
        if batch_size > 0:
            success = run(0, batch_size)
            if success != 0:
                raise ComputationError(f"PFAPACK returned error code {success}")
    else:
        # ctypes releases the GIL during the C call and every chunk is an
        # independent contiguous view, so threads scale near-linearly.
        bounds = np.linspace(0, batch_size, workers + 1, dtype=np.intp)
        with ThreadPoolExecutor(max_workers=workers) as executor:
            results = executor.map(run, bounds[:-1], bounds[1:])
        for success in results:
            if success != 0:
                raise ComputationError(f"PFAPACK returned error code {success}")

    # For 2D input batch_shape is (), so this returns a scalar.
    return pfaffians.reshape(batch_shape)[()]
