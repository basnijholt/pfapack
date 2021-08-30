import numpy as np
import pfapack.pfaffian as pf


def int_rand_mat(n):
    A = np.matrix(list(map(np.vectorize(round), np.random.rand(n, n) * 10)))
    return A - A.T


def float_rand_mat(n):
    A = np.random.rand(n, n) * 10
    return A - A.T


def complex_rand_mat(n):
    A = np.random.rand(n, n) * 10 + np.random.rand(n, n) * 10j
    return A - A.T


def int_rand_mat2(n):
    A = np.matrix(list(map(np.around, np.random.rand(n, n) * 10)))
    return A - A.T


def test_integer_matrices():
    for i in range(2, 12, 2):
        A = int_rand_mat(i)
        H = pf.pfaffian(A, method="H")
        P = pf.pfaffian(A, method="P")
        result_H = np.abs(H ** 2 - np.linalg.det(A))
        result_P = np.abs(P ** 2 - np.linalg.det(A))
        eps = 10 ** (-i)
        assert result_H < eps
        assert result_P < eps


def test_real_matrices():
    for i in range(2, 12, 2):
        A = float_rand_mat(i)
        H = pf.pfaffian(A, method="H")
        P = pf.pfaffian(A, method="P")
        result_H = np.abs(H ** 2 - np.linalg.det(A))
        result_P = np.abs(P ** 2 - np.linalg.det(A))
        eps = 10 ** (-i)
        assert result_H < eps
        assert result_P < eps


def test_complex_matrices():
    for i in range(2, 12, 2):
        A = complex_rand_mat(i)
        H = pf.pfaffian(A, method="H")
        P = pf.pfaffian(A, method="P")
        result_H = np.abs(H ** 2 - np.linalg.det(A))
        result_P = np.abs(P ** 2 - np.linalg.det(A))
        eps = 10 ** (-i)
        assert result_H < eps
        assert result_P < eps
