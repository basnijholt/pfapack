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


print("Testing integer matrices...")
for i in np.arange(2, 12, 2):
    A = int_rand_mat(6)
    H = pf.pfaffian(A, method="H")
    P = pf.pfaffian(A, method="P")
    print("Method H: ", H ** 2 - np.linalg.det(A))
    print("Method P: ", P ** 2 - np.linalg.det(A))

print("Testing real matrices...")
for i in np.arange(2, 12, 2):
    A = float_rand_mat(6)
    H = pf.pfaffian(A, method="H")
    P = pf.pfaffian(A, method="P")
    print("Method H: ", H ** 2 - np.linalg.det(A))
    print("Method P: ", P ** 2 - np.linalg.det(A))


print("Testing complex matrices...")
for i in np.arange(2, 12, 2):
    A = complex_rand_mat(6)
    H = pf.pfaffian(A, method="H")
    P = pf.pfaffian(A, method="P")
    print("Method H: ", H ** 2 - np.linalg.det(A))
    print("Method P: ", P ** 2 - np.linalg.det(A))
