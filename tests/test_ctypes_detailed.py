from __future__ import annotations

import numpy as np
import numpy.matlib
import pytest

from pfapack.ctypes import pfaffian as cpf
from pfapack.exceptions import InvalidDimensionError, InvalidParameterError

try:
    from pfapack.ctypes import pfaffian as cpf

    with_ctypes = True
except OSError:
    with_ctypes = False


@pytest.mark.skipif(not with_ctypes, reason="the libs might not be installed")
def test_ctypes_real_different_sizes():
    """Test real matrices of different sizes."""
    for n in [2, 4, 8, 16, 32]:
        A = numpy.matlib.rand(n, n)
        A = A - A.T  # make skew-symmetric

        # Test different methods
        pf_p = cpf(A, uplo="L", method="P")
        pf_h = cpf(A, uplo="L", method="H")

        # Test overflow avoidance
        pf_p_safe = cpf(A, uplo="L", method="P", avoid_overflow=True)
        pf_h_safe = cpf(A, uplo="L", method="H", avoid_overflow=True)

        # All methods should give the same result
        np.testing.assert_almost_equal(pf_p / pf_h, 1, decimal=10)
        np.testing.assert_almost_equal(pf_p / pf_p_safe, 1, decimal=10)
        np.testing.assert_almost_equal(pf_h / pf_h_safe, 1, decimal=10)

        # Test that Pf² = det(A)
        det_a = np.linalg.det(A)
        np.testing.assert_almost_equal(pf_p * pf_p / det_a, 1, decimal=10)


@pytest.mark.skipif(not with_ctypes, reason="the libs might not be installed")
def test_ctypes_complex_different_sizes():
    """Test complex matrices of different sizes."""
    for n in [2, 4, 8, 16, 32]:
        A = numpy.matlib.rand(n, n) + 1.0j * numpy.matlib.rand(n, n)
        A = A - A.T  # make skew-symmetric

        # Test different methods
        pf_p = cpf(A, uplo="L", method="P")
        pf_h = cpf(A, uplo="L", method="H")

        # Test overflow avoidance
        pf_p_safe = cpf(A, uplo="L", method="P", avoid_overflow=True)
        pf_h_safe = cpf(A, uplo="L", method="H", avoid_overflow=True)

        # All methods should give the same result
        np.testing.assert_almost_equal(pf_p / pf_h, 1, decimal=10)
        np.testing.assert_almost_equal(pf_p / pf_p_safe, 1, decimal=10)
        np.testing.assert_almost_equal(pf_h / pf_h_safe, 1, decimal=10)

        # Test that Pf² = det(A)
        det_a = np.linalg.det(A)
        np.testing.assert_almost_equal(pf_p * pf_p / det_a, 1, decimal=10)


@pytest.mark.skipif(not with_ctypes, reason="the libs might not be installed")
def test_ctypes_uplo():
    """Test that upper and lower triangular options give same results."""
    n = 8
    A = numpy.matlib.rand(n, n)
    A = A - A.T

    pf_upper = cpf(A, uplo="U", method="P")
    pf_lower = cpf(A, uplo="L", method="P")

    np.testing.assert_almost_equal(pf_upper / pf_lower, 1, decimal=10)


def make_skew_matrix(n: int, complex: bool = False) -> np.ndarray:
    """Create a random skew-symmetric matrix."""
    A = np.random.rand(n, n)
    if complex:
        A = A + 1j * np.random.rand(n, n)
    return A - A.T


@pytest.mark.skipif(not with_ctypes, reason="the libs might not be installed")
def test_ctypes_errors():
    """Test error conditions."""
    # Test non-square matrix
    A = np.random.rand(3, 4)
    with pytest.raises(InvalidDimensionError, match="must be square"):
        cpf(A)

    # Test odd-dimensional matrix
    A = make_skew_matrix(3)
    with pytest.raises(InvalidDimensionError, match="must be even"):
        cpf(A)

    # Test invalid uplo parameter
    A = make_skew_matrix(4)
    with pytest.raises(InvalidParameterError, match="'uplo' must be"):
        cpf(A, uplo="X")

    # Test invalid method parameter
    with pytest.raises(InvalidParameterError, match="'method' must be"):
        cpf(A, method="X")
