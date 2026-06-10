import numpy as np
import numpy.linalg
import numpy.random
import pytest

from pfapack import pfaffian as pf  # noqa isort:skip

try:
    from pfapack.ctypes import pfaffian as cpfaffian

    with_ctypes = True
except OSError:
    with_ctypes = False


EPS = 1e-11


def test_pfaffian():
    # Compare the output of the different Pfaffian routines
    # and compare to the determinant

    # first real matrices
    A = numpy.random.rand(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method="H")
    pfa3 = pf.pfaffian_schur(A)
    print(pfa1, pfa2, pfa3)
    deta = numpy.linalg.det(A)

    assert abs((pfa1 - pfa2) / pfa1) < EPS
    assert abs((pfa1 - pfa3) / pfa1) < EPS
    assert abs((pfa1**2 - deta) / deta) < EPS

    # then complex matrices
    A = numpy.random.rand(100, 100) + 1.0j * numpy.random.rand(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method="H")

    deta = numpy.linalg.det(A)

    assert abs((pfa1 - pfa2) / pfa1) < EPS
    assert abs((pfa1**2 - deta) / deta) < EPS


def test_decompositions():
    # Test the LTL^T and Householder decompositions

    # first real matrices
    A = numpy.random.rand(100, 100)
    A = A - A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P * A * P.T - L * T * L.T) / numpy.linalg.norm(A) < EPS

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A - Q * T * Q.T) / numpy.linalg.norm(A) < EPS

    # then complex matrices
    A = numpy.random.rand(100, 100) + 1.0j * numpy.random.rand(100, 100)
    A = A - A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P * A * P.T - L * T * L.T) / numpy.linalg.norm(A) < EPS

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A - Q * T * Q.T) / numpy.linalg.norm(A) < EPS


def test_skew_symmetry_check_negative_deviation():
    # Regression test: the old check `abs((A + A.T).max()) < 1e-14` only
    # looked at the largest entry of A + A.T, so a matrix whose deviation
    # from skew-symmetry is purely *negative* slipped through silently.
    A = np.array([[0.0, -1.0], [-1.0, 0.0]])  # symmetric, (A + A.T).max() == 0

    for func in (
        pf.pfaffian,
        pf.pfaffian_LTL,
        pf.pfaffian_householder,
        pf.pfaffian_schur,
        pf.skew_LTL,
        pf.skew_tridiagonalize,
    ):
        with pytest.raises(AssertionError):
            func(A)


def test_skew_symmetry_check_complex():
    # The old check compared complex numbers lexicographically (real part
    # first), which is not a magnitude check; this matrix passed it.
    A = np.array([[0.0, -1.0 + 1.0j], [-1.0 + 1.0j, 0.0]])

    for func in (
        pf.pfaffian,
        pf.pfaffian_LTL,
        pf.pfaffian_householder,
        pf.skew_LTL,
        pf.skew_tridiagonalize,
    ):
        with pytest.raises(AssertionError):
            func(A)


def test_skew_symmetry_check_large_scale():
    # A skew-symmetric matrix with entries ~1e8 carries roundoff of order
    # 1e-8 in A + A.T; an absolute tolerance would spuriously reject it.
    rng = numpy.random.default_rng(0)
    A = rng.random((10, 10))
    A = (A - A.T) * 1e8
    # simulate floating-point roundoff at the matrix scale
    A += rng.standard_normal((10, 10)) * 1e-8

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method="H")
    pfa3 = pf.pfaffian_schur(A)

    assert abs((pfa1 - pfa2) / pfa1) < EPS
    assert abs((pfa1 - pfa3) / pfa1) < EPS


def test_skew_symmetry_check_tiny_scale():
    # A tiny-scaled skew-symmetric matrix must pass the check and yield
    # the correct (tiny) Pfaffian.
    rng = numpy.random.default_rng(0)
    A = rng.random((8, 8))
    A = (A - A.T) * 1e-12

    pfa1 = pf.pfaffian_LTL(A)
    pfa2 = pf.pfaffian_householder(A)

    assert pfa1 != 0
    assert abs((pfa1 - pfa2) / pfa1) < EPS


@pytest.mark.skipif(not with_ctypes, reason="the libs might not be installed")
def test_ctypes():
    for method in ("P", "H"):
        # first real matrices
        A = numpy.random.rand(100, 100)
        A = A - A.T
        pf_a = cpfaffian(A, uplo="L", method=method)
        pf_a2 = cpfaffian(A, uplo="L", avoid_overflow=True, method=method)

        np.testing.assert_almost_equal(pf_a / pf_a2, 1)
        np.testing.assert_almost_equal(pf_a / pf.pfaffian(A), 1)

        # then complex matrices
        A = numpy.random.rand(100, 100) + 1.0j * numpy.random.rand(100, 100)
        A = A - A.T
        pf_a = cpfaffian(A, uplo="L", method=method)
        pf_a2 = cpfaffian(A, uplo="L", avoid_overflow=True, method=method)

        np.testing.assert_almost_equal(pf_a / pf_a2, 1)
        np.testing.assert_almost_equal(pf_a / pf.pfaffian(A), 1)
