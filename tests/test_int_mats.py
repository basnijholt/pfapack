import numpy as np
import pytest

import pfapack.pfaffian as pf


def int_rand_mat(n):
    A = np.array(list(map(np.vectorize(round), np.random.rand(n, n))))
    return A - A.T


def float_rand_mat(n):
    A = np.random.rand(n, n)
    return A - A.T


def complex_rand_mat(n):
    A = np.random.rand(n, n) + 1j * np.random.rand(n, n)
    return A - A.T


def int_rand_mat2(n):
    A = np.array(list(map(np.around, np.random.rand(n, n))))
    return A - A.T


@pytest.mark.parametrize(
    "f_matrix",
    [int_rand_mat, float_rand_mat, complex_rand_mat, int_rand_mat2],
)
def test_matrices(f_matrix):
    for i in range(40, 200, 20):
        A = f_matrix(i)
        H = pf.pfaffian(A, method="H")
        P = pf.pfaffian(A, method="P")
        da = np.linalg.det(A)
        result_H = np.log10(H**2) - np.log10(da)
        result_P = np.log10(P**2) - np.log10(da)
        assert result_H < 1e-11
        assert result_P < 1e-11
