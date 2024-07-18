"""Cython version of pfapack (just Parlett-Reid algorithm at the moment)."""
cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport fabs, sqrt
from libc.complex cimport cabs


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] outer_real(double[:] a, double[:] b):
    cdef Py_ssize_t dim_a = a.size
    cdef Py_ssize_t dim_b = b.size
    cdef np.ndarray out = np.empty((dim_a, dim_b), dtype=np.float64)
    cdef double[:, :] out_view = out

    cdef int i, j
    for i in range(dim_a):
        for j in range(dim_b):
            out_view[i, j] = a[i] * b[j]

    return out_view

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double complex[:, :] outer_complex(double complex[:] a, double complex[:] b):
    cdef Py_ssize_t dim_a = a.size
    cdef Py_ssize_t dim_b = b.size
    cdef np.ndarray[np.complex128_t, ndim=2] out = np.empty((dim_a, dim_b), dtype=np.complex128)
    cdef double complex[:, :] out_view = out

    cdef int i, j
    for i in range(dim_a):
        for j in range(dim_b):
            out_view[i, j] = a[i] * b[j]

    return out_view

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double pfaffian_real(double[:, :] A_view, bint overwrite_a=False):
    cdef:
        Py_ssize_t n = A_view.shape[0]
        Py_ssize_t k, kp, i, j
        double max_value, pfaffian_val = 1.0, temp
        double[:] tau

    if n % 2 == 1:
        return 0.0

    if not overwrite_a:
        A_copy = np.array(A_view, copy=True)
        A_view = A_copy

    for k in range(0, n - 1, 2):
        max_index = k + 1
        max_value = fabs(A_view[k+1, k])
        for i in range(k + 2, n):
            if fabs(A_view[i, k]) > max_value:
                max_value = fabs(A_view[i, k])
                max_index = i
        kp = max_index

        if kp != k + 1:
            for j in range(k, n):
                temp = A_view[k + 1, j]
                A_view[k + 1, j] = A_view[kp, j]
                A_view[kp, j] = temp
            for i in range(k, n):
                temp = A_view[i, k + 1]
                A_view[i, k + 1] = A_view[i, kp]
                A_view[i, kp] = temp
            pfaffian_val *= -1

        if A_view[k, k + 1] != 0.0:
            pfaffian_val *= A_view[k, k + 1]
            if k + 2 < n:
                tau = np.empty(n - (k + 2), dtype=np.float64)
                for i in range(n - (k + 2)):
                    tau[i] = A_view[k, k + 2 + i] / A_view[k, k + 1]
                addition = outer_real(tau, A_view[k + 2:, k + 1])
                subtraction = outer_real(A_view[k + 2:, k + 1], tau)

                # Manual update without using inplace operations
                for i in range(n - (k + 2)):
                    for j in range(n - (k + 2)):
                        A_view[k + 2 + i, k + 2 + j] += addition[i, j]
                        A_view[k + 2 + i, k + 2 + j] -= subtraction[i, j]

    return pfaffian_val

cpdef pfaffian(np.ndarray A, bint overwrite_a=False):
    """
    Compute the Pfaffian of a real or complex skew-symmetric
    matrix A (A=-A^T). If overwrite_a=True, the matrix A
    is overwritten in the process.
    
    Depending on the data type of A, it calls either pfaffian_real or
    pfaffian_complex functions.
    """
    cdef int n = A.shape[0]
    
    # Check if matrix is square
    assert n == A.shape[1] and n > 0
    # Check if it's skew-symmetric
    assert abs((A + A.T).max()) < 1e-14
    
    if np.iscomplexobj(A):
        return pfaffian_complex(A, overwrite_a)
    else:
        return pfaffian_real(A, overwrite_a)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef double complex pfaffian_complex(double complex[:, :] A_view, bint overwrite_a=False):
    cdef:
        Py_ssize_t n = A_view.shape[0]
        Py_ssize_t k, kp, i, j
        double max_value = 1.0
        double complex pfaffian_val = 1.0 + 0.0j
        double complex temp
        double complex[:] tau

    if n % 2 == 1:
        return 0.0 + 0.0j

    if not overwrite_a:
        A_copy = np.array(A_view, copy=True)
        A_view = A_copy

    for k in range(0, n - 1, 2):
        max_index = k + 1
        max_value = cabs(A_view[k+1, k])
        for i in range(k + 2, n):
            if cabs(A_view[i, k]) > max_value:
                max_value = cabs(A_view[i, k])
                max_index = i
        kp = max_index

        if kp != k + 1:
            for j in range(k, n):
                temp = A_view[k + 1, j]
                A_view[k + 1, j] = A_view[kp, j]
                A_view[kp, j] = temp
            for i in range(k, n):
                temp = A_view[i, k + 1]
                A_view[i, k + 1] = A_view[i, kp]
                A_view[i, kp] = temp
            pfaffian_val *= -1

        if A_view[k, k + 1] != 0.0 + 0.0j:
            pfaffian_val *= A_view[k, k + 1]
            if k + 2 < n:
                tau = np.empty(n - (k + 2), dtype=np.complex128)
                for i in range(n - (k + 2)):
                    tau[i] = A_view[k, k + 2 + i] / A_view[k, k + 1]
                addition = outer_complex(tau, A_view[k + 2:, k + 1])
                subtraction = outer_complex(A_view[k + 2:, k + 1], tau)

                # Manual update without using inplace operations
                for i in range(n - (k + 2)):
                    for j in range(n - (k + 2)):
                        A_view[k + 2 + i, k + 2 + j] += addition[i, j]
                        A_view[k + 2 + i, k + 2 + j] -= subtraction[i, j]

    return pfaffian_val

def test_pfaffian():
    print("")
    print("test_pfaffian:")
    print("  Compare the output of Pfaffian routines and determinant.")
    #
    #  Real matrices
    #
    A = np.random.randn(100, 100)
    A = A - A.T

    pfa1 = pfaffian(A)
    deta = np.linalg.det(A)

    print("")
    print("  Real matrix:")
    print("    pfaffian(A) =           ", pfa1)
    print("")
    print("    pfaffian(A)^2 =         ", pfa1**2)
    print("    np.linalg.det(A) =   ", deta)

    assert abs((pfa1**2 - deta) / deta) < 1e-13
    #
    #  Complex matrices
    #
    A = np.random.randn(100, 100) + 1.0j * np.random.randn(100, 100)
    A = A - A.T

    pfa1 = pfaffian(A)
    deta = np.linalg.det(A)

    print("")
    print("  Complex matrix:")
    print("    pfaffian(A) =            ", pfa1)
    print("")
    print("    pfaffian(A)^2 =          ", pfa1**2)
    print("    np.linalg.det(A) =    ", deta)

    assert abs((pfa1**2 - deta) / deta) < 1e-13

    return


def timestamp():
    import time

    t = time.time()
    print(time.ctime(t))

    return


def toms923_test():
    import platform

    print("")
    print("toms923_test():")
    print("  Python version: %s" % (platform.python_version()))
    print("  Test toms923().")
    test_pfaffian()
    print("")
    print("toms923_test():")
    print("  Normal end of execution.")
    return


if __name__ == "__main__":
    timestamp()
    toms923_test()
    timestamp()
