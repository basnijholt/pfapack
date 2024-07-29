import numpy as np
import numpy.linalg
import pytest

import sys  # isort:skip

sys.path.append("..")

from pfapack import pfaffian as pf  # noqa isort:skip

from pfapack.ctypes import pfaffian as cpfaffian
from pfapack.ctypes import pfaffian_batched as cpfaffian_batched


EPS = 1e-11


def test_pfaffian():
    # Compare the output of the different Pfaffian routines
    # and compare to the determinant

    # first real matrices
    A = np.random.randn(100, 100)
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
    A = np.random.randn(100, 100) + 1.0j * np.random.randn(100, 100)
    A = A - A.T

    pfa1 = pf.pfaffian(A)
    pfa2 = pf.pfaffian(A, method="H")

    deta = numpy.linalg.det(A)

    assert abs((pfa1 - pfa2) / pfa1) < EPS
    assert abs((pfa1 ** 2 - deta) / deta) < EPS


def test_decompositions():
    # Test the LTL^T and Householder decompositions

    # first real matrices
    A = np.random.randn(100, 100)
    A = A - A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P * A * P.T - L * T * L.T) / numpy.linalg.norm(A) < EPS

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A - Q * T * Q.T) / numpy.linalg.norm(A) < EPS

    # then complex matrices
    A = np.random.randn(100, 100) + 1.0j * np.random.randn(100, 100)
    A = A - A.T

    T, L, P = pf.skew_LTL(A)

    assert numpy.linalg.norm(P * A * P.T - L * T * L.T) / numpy.linalg.norm(A) < EPS

    T, Q = pf.skew_tridiagonalize(A)

    assert numpy.linalg.norm(A - Q * T * Q.T) / numpy.linalg.norm(A) < EPS


def test_ctypes():
    for method in ("P", "H"):
        # first real matrices
        A = np.random.randn(100, 100)
        A = A - A.T
        pf_a = cpfaffian(A, uplo="L", method=method)
        pf_a2 = cpfaffian(A, uplo="L", avoid_overflow=True, method=method)

        np.testing.assert_almost_equal(pf_a / pf_a2, 1)
        np.testing.assert_almost_equal(pf_a / pf.pfaffian(A), 1)

        # then complex matrices
        A = np.random.randn(100, 100) + 1.0j * np.random.randn(100, 100)
        A = A - A.T
        pf_a = cpfaffian(A, uplo="L", method=method)
        pf_a2 = cpfaffian(A, uplo="L", avoid_overflow=True, method=method)

        np.testing.assert_almost_equal(pf_a / pf_a2, 1)
        np.testing.assert_almost_equal(pf_a / pf.pfaffian(A), 1)

def test_batched_vs_individual_float64():
    batch_size = 20
    dtype = np.float64
    for matrix_size in [4,6,8,10,12,14,16]:
        # Generate a batch of random skew-symmetric matrices
        batch = np.random.randn(batch_size, matrix_size, matrix_size).astype(dtype)
        batch = batch - np.transpose(batch, (0, 2, 1))

        # Calculate Pfaffians using batched method
        pfaffians_batched = cpfaffian_batched(batch)

        # Calculate Pfaffians individually
        pfaffians_individual = np.array([cpfaffian(matrix) for matrix in batch])

        # Compare results
        np.testing.assert_allclose(pfaffians_batched, pfaffians_individual, rtol=EPS, atol=EPS)

def test_batched_vs_individual_complex128():
    batch_size = 20
    dtype = np.complex128
    for matrix_size in [4, 6, 8, 10, 12, 14, 16]:
        # Generate a batch of random skew-Hermitian matrices
        real_part = np.random.randn(batch_size, matrix_size, matrix_size)
        imag_part = np.random.randn(batch_size, matrix_size, matrix_size)
        batch = (real_part + 1j * imag_part).astype(dtype)
        # Make the matrices skew-Hermitian
        batch = batch - np.conjugate(np.transpose(batch, (0, 2, 1)))
        # Calculate Pfaffians using batched method
        pfaffians_batched = cpfaffian_batched(batch)
        # Calculate Pfaffians individually
        pfaffians_individual = np.array([cpfaffian(matrix) for matrix in batch])
        # Compare results
        assert np.allclose(pfaffians_batched, pfaffians_individual, rtol=EPS, atol=EPS)

def test_known_values():
    # Test with known values for both real and complex matrices
    # Real skew-symmetric matrix
    real_matrix = np.array([
        [ 0,  1, -2,  3],
        [-1,  0,  4, -5],
        [ 2, -4,  0,  6],
        [-3,  5, -6,  0]
    ], dtype=np.float64)
    expected_real_pfaffian = pf.pfaffian(real_matrix) 

    # Complex skew-Hermitian matrix
    complex_matrix = np.array([
        [ 0,    1-1j, -2+2j,  3-3j],
        [-1+1j,  0,    4-4j, -5+5j],
        [ 2-2j, -4+4j,  0,    6-6j],
        [-3+3j,  5-5j, -6+6j,  0   ]
    ], dtype=np.complex128)
    expected_complex_pfaffian = pf.pfaffian(complex_matrix) 

    assert np.allclose(cpfaffian(real_matrix), expected_real_pfaffian, rtol=EPS, atol=EPS)
    assert np.allclose(cpfaffian(complex_matrix), expected_complex_pfaffian, rtol=EPS, atol=EPS)

    # Test batched version with known values
    batch_real = np.array([real_matrix, real_matrix])
    batch_complex = np.array([complex_matrix, complex_matrix])

    assert np.allclose(cpfaffian_batched(batch_real), 
                    np.array([expected_real_pfaffian, expected_real_pfaffian]), 
                    rtol=EPS, atol=EPS)
    assert np.allclose(cpfaffian_batched(batch_complex), 
                    np.array([expected_complex_pfaffian, expected_complex_pfaffian]), 
                    rtol=EPS, atol=EPS)

