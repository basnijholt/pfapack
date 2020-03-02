import numpy.linalg
import numpy.matlib

import sys  # isort:skip

sys.path.append("..")

from pfapack import pfaffian as pf  # noqa isort:skip

EPS = 1e-12


def test_pfaffian():
    # Compare the output of the different Pfaffian routines
    # and compare to the determinant

    # first real matrices
    A = numpy.matlib.rand(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method="H")
    pfa3 = pf.pfaffian_schur(A)
    print(pfa1, pfa2, pfa3)
    deta = numpy.linalg.det(A)

    assert abs((pfa1 - pfa2) / pfa1) < EPS
    assert abs((pfa1 - pfa3) / pfa1) < EPS
    assert abs((pfa1 ** 2 - deta) / deta) < EPS

    # then complex matrices
    A = numpy.matlib.rand(100, 100) + 1.0j * numpy.matlib.rand(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method="H")

    deta = numpy.linalg.det(A)

    assert abs((pfa1 - pfa2) / pfa1) < EPS
    assert abs((pfa1 ** 2 - deta) / deta) < EPS


def test_decompositions():
    # Test the LTL^T and Householder decompositions

    # first real matrices
    A = numpy.matlib.rand(100, 100)
    A = A - A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P * A * P.T - L * T * L.T) / numpy.linalg.norm(A) < EPS

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A - Q * T * Q.T) / numpy.linalg.norm(A) < EPS

    # then complex matrices
    A = numpy.matlib.rand(100, 100) + 1.0j * numpy.matlib.rand(100, 100)
    A = A - A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P * A * P.T - L * T * L.T) / numpy.linalg.norm(A) < EPS

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A - Q * T * Q.T) / numpy.linalg.norm(A) < EPS
