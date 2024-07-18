"""A package for computing Pfaffians"""
cimport cython

import cmath
import math

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

from libc.math cimport sqrt, atan2, exp, fabs
from libc.complex cimport csqrt, cexp, conj
cimport numpy as cnp
cimport cython


def get_column_slice_as_1d_real(cnp.ndarray[cnp.float64_t, ndim=2] A, int start_row, int col):
    """
    Extract a column slice as a 1D array from the given real matrix.
    """
    cdef:
        int n = A.shape[0] - start_row
        cnp.ndarray[cnp.float64_t, ndim=1] slice_1d = np.empty(n, dtype=A.dtype)
        int i
    for i in range(n):
        slice_1d[i] = A[start_row + i, col]
    return slice_1d

def get_column_slice_as_1d_complex(cnp.ndarray[cnp.complex128_t, ndim=2] A, int start_row, int col):
    """
    Extract a column slice as a 1D array from the given complex matrix.
    """
    cdef:
        int n = A.shape[0] - start_row
        cnp.ndarray[cnp.complex128_t, ndim=1] slice_1d = np.empty(n, dtype=A.dtype)
        int i

    for i in range(n):
        slice_1d[i] = A[start_row + i, col]
    return slice_1d

@cython.boundscheck(False)
@cython.wraparound(False)
def householder_real(cnp.ndarray[cnp.float64_t, ndim=1] x):
    """(v, tau, alpha) = householder_real(x)

    Compute a Householder transformation such that
    (1-tau v v^T) x = alpha e_1
    where x and v a real vectors, tau is 0 or 2, and
    alpha a real number (e_1 is the first unit vector)
    """

    cdef:
        cnp.float64_t sigma = 0
        cnp.float64_t norm_x = 0
        cnp.float64_t alpha = 0
        int i
        cnp.ndarray[cnp.float64_t, ndim=1] v

    # Compute sigma, the sum of squares of x[1:]
    for i in range(1, x.shape[0]):
        sigma += x[i] * x[i]

    if sigma == 0:
        v = np.zeros(x.shape[0], dtype=np.float64)
        return (v, 0, x[0])
    else:
        norm_x = sqrt(x[0] ** 2 + sigma)
        v = x.copy()

        # Adjust the sign of the first element and compute alpha
        if x[0] <= 0:
            v[0] -= norm_x
            alpha = norm_x
        else:
            v[0] += norm_x
            alpha = -norm_x

        # Normalize v
        norm_v = sqrt(np.dot(v, v))
        for i in range(v.shape[0]):
            v[i] /= norm_v

        return (v, 2, alpha)

@cython.boundscheck(False)
def householder_complex(cnp.ndarray[cnp.complex128_t, ndim=1] x):
    """(v, tau, alpha) = householder_complex(x)

    Compute a Householder transformation such that
    (1-tau v v^T) x = alpha e_1
    where x and v are complex vectors, tau is 0 or 2, and
    alpha is a complex number (e_1 is the first unit vector)
    """
    assert x.shape[0] > 0

    cdef cnp.complex128_t sigma = 0 + 0j
    cdef cnp.complex128_t norm_x = 0 + 0j
    cdef cnp.complex128_t phase = 0 + 0j
    cdef int i
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] v

    # Compute sigma as the dot product of x[1:] with its conjugate
    for i in range(1, len(x)):
        sigma += x[i].conjugate() * x[i]

    if sigma == 0:
        return (np.zeros(x.shape[0], dtype=np.complex128), 0, x[0])
    else:
        norm_x = csqrt(x[0].conjugate() * x[0] + sigma)

        v = x.copy()

        phase = cexp(1j * atan2(x[0].imag, x[0].real))

        v[0] += phase * norm_x

        # Normalize v
        v /= np.linalg.norm(v)

    return (v, 2, -phase * norm_x)

@cython.boundscheck(False)
def skew_tridiagonalize_real(cnp.ndarray[cnp.float64_t, ndim=2] A, bint overwrite_a=False, bint calc_q=True):
    """
    Transform a skew-symmetric real matrix into tridiagonal form with orthogonal transformation.
    """

    cdef:
        int i, n = A.shape[0]
        cnp.ndarray[cnp.float64_t, ndim=1] v
        cnp.float64_t tau, alpha, w_elem
        cnp.ndarray[cnp.float64_t, ndim=1] w, y
        cnp.ndarray[cnp.float64_t, ndim=2] Q

    # Basic checks
    assert n == A.shape[1] and n > 0, "Matrix must be square and non-empty"
    assert abs((A + A.T).max()) < 1e-14, "Matrix must be skew-symmetric"

    if not overwrite_a:
        A = A.copy()

    if calc_q:
        Q = np.eye(n, dtype=A.dtype)

    for i in range(n - 2):
        v = get_column_slice_as_1d_real(A, i + 1, i)
        v, tau, alpha = householder_real(v)
        A[i + 1, i] = alpha
        A[i, i + 1] = -alpha
        A[i + 2 :, i] = 0
        A[i, i + 2 :] = 0

        w = tau * np.dot(A[i + 1 :, i + 1 :], v)
        A[i + 1 :, i + 1 :] += np.outer(v, w) - np.outer(w, v)

        if calc_q:
            y = tau * np.dot(Q[:, i + 1 :], v)
            Q[:, i + 1 :] -= np.outer(y, v)

    if calc_q:
        return (np.asmatrix(A), np.asmatrix(Q))
    else:
        return np.asmatrix(A)


@cython.boundscheck(False)
def skew_tridiagonalize_complex(cnp.ndarray[cnp.complex128_t, ndim=2] A, bint overwrite_a=False, bint calc_q=True):
    """
    Transform a skew-symmetric complex matrix into tridiagonal form with unitary transformation.
    """

    cdef:
        int i, n = A.shape[0]
        cnp.ndarray[cnp.complex128_t, ndim=1] v
        cnp.complex128_t tau, alpha, w_elem
        cnp.ndarray[cnp.complex128_t, ndim=1] w, y
        cnp.ndarray[cnp.complex128_t, ndim=2] Q

    # Basic checks
    assert n == A.shape[1] and n > 0, "Matrix must be square and non-empty"
    assert abs((A + A.T).max()) < 1e-14, "Matrix must be skew-symmetric"

    if not overwrite_a:
        A = A.copy()

    if calc_q:
        Q = np.eye(n, dtype=A.dtype)

    for i in range(n - 2):
        v = get_column_slice_as_1d_complex(A, i + 1, i)
        v, tau, alpha = householder_complex(v)
        A[i + 1, i] = alpha
        A[i, i + 1] = -alpha
        A[i + 2 :, i] = 0
        A[i, i + 2 :] = 0

        # Ensure v is contiguous
        v_conj = np.ascontiguousarray(v.conj())
        A_sub = np.ascontiguousarray(A[i+1:, i+1:])

        # Update the matrix block A(i+1:N, i+1:N)
        w = tau * np.dot(A_sub, v_conj)
        A[i+1:, i+1:] += np.outer(v, w) - np.outer(w, v)

        if calc_q:
            y = tau * np.dot(Q[:, i + 1 :], v)
            Q[:, i + 1 :] -= np.outer(y, np.conj(v))

    if calc_q:
        return (np.asmatrix(A), np.asmatrix(Q))
    else:
        return np.asmatrix(A)


@cython.boundscheck(False)
def skew_tridiagonalize(cnp.ndarray A, bint overwrite_a=False, bint calc_q=True):
    """
    Driver function to transform a skew-symmetric matrix into tridiagonal form.
    Selects the appropriate function based on the dtype of the matrix.
    """

    if np.issubdtype(A.dtype, np.complexfloating):
        return skew_tridiagonalize_complex(A, overwrite_a, calc_q)
    elif np.issubdtype(A.dtype, np.floating):
        return skew_tridiagonalize_real(A, overwrite_a, calc_q)
    else:
        raise TypeError("skew_tridiagonalize() can only work on numeric input")

def skew_LTL(A, overwrite_a=False, calc_L=True, calc_P=True):
    """T, L, P = skew_LTL(A, overwrite_a, calc_q=True)

    Bring a real or complex skew-symmetric matrix (A=-A^T) into
    tridiagonal form T (with zero diagonal) with a lower unit
    triangular matrix L such that
    P A P^T= L T L^T

    A is overwritten if overwrite_a=True (default: False),
    L and P only calculated if calc_L=True or calc_P=True,
    respectively (default: True).
    """

    # Check if matrix is square
    assert A.shape[0] == A.shape[1] > 0
    # Check if it's skew-symmetric
    assert abs((A + A.T).max()) < 1e-14

    n = A.shape[0]
    A = np.asarray(A)  # the slice views work only properly for arrays

    if not overwrite_a:
        A = A.copy()

    if calc_L:
        L = np.eye(n, dtype=A.dtype)

    if calc_P:
        Pv = np.arange(n)

    for k in range(n - 2):
        # First, find the largest entry in A[k+1:,k] and
        # permute it to A[k+1,k]
        kp = k + 1 + np.abs(A[k + 1 :, k]).argmax()

        # Check if we need to pivot
        if kp != k + 1:
            # interchange rows k+1 and kp
            temp = A[k + 1, k:].copy()
            A[k + 1, k:] = A[kp, k:]
            A[kp, k:] = temp

            # Then interchange columns k+1 and kp
            temp = A[k:, k + 1].copy()
            A[k:, k + 1] = A[k:, kp]
            A[k:, kp] = temp

            if calc_L:
                # permute L accordingly
                temp = L[k + 1, 1 : k + 1].copy()
                L[k + 1, 1 : k + 1] = L[kp, 1 : k + 1]
                L[kp, 1 : k + 1] = temp

            if calc_P:
                # accumulate the permutation matrix
                temp = Pv[k + 1]
                Pv[k + 1] = Pv[kp]
                Pv[kp] = temp

        # Now form the Gauss vector
        if A[k + 1, k] != 0.0:
            tau = A[k + 2 :, k].copy()
            tau /= A[k + 1, k]

            # clear eliminated row and column
            A[k + 2 :, k] = 0.0
            A[k, k + 2 :] = 0.0

            # Update the matrix block A(k+2:,k+2)
            A[k + 2 :, k + 2 :] += np.outer(tau, A[k + 2 :, k + 1])
            A[k + 2 :, k + 2 :] -= np.outer(A[k + 2 :, k + 1], tau)

            if calc_L:
                L[k + 2 :, k + 1] = tau

    if calc_P:
        # form the permutation matrix as a sparse matrix
        P = sp.csr_matrix((np.ones(n), (np.arange(n), Pv)))

    if calc_L:
        if calc_P:
            return (np.asmatrix(A), np.asmatrix(L), P)
        else:
            return (np.asmatrix(A), np.asmatrix(L))
    else:
        if calc_P:
            return (np.asmatrix(A), P)
        else:
            return np.asmatrix(A)

@cython.boundscheck(False)
@cython.wraparound(False)
def pfaffian(cnp.ndarray A, bint overwrite_a=False, str method='P'):
    assert A.shape[0] == A.shape[1] > 0, "Matrix must be square and non-empty"
    assert fabs((A + A.T).max()) < 1e-14, "Matrix must be skew-symmetric"
    assert method == 'P' or method == 'H', "Method must be 'P' or 'H'"

    if method == 'P':
        return pfaffian_LTL(A, overwrite_a)
    else:
        return pfaffian_householder(A, overwrite_a)

def pfaffian_LTL(A, overwrite_a=False):
    """pfaffian_LTL(A, overwrite_a=False)

    Compute the Pfaffian of a real or complex skew-symmetric
    matrix A (A=-A^T). If overwrite_a=True, the matrix A
    is overwritten in the process. This function uses
    the Parlett-Reid algorithm.
    """
    # Check if matrix is square
    assert A.shape[0] == A.shape[1] > 0
    # Check if it's skew-symmetric
    assert abs((A + A.T).max()) < 1e-14

    n, m = A.shape
    # type check to fix problems with integer numbers
    dtype = type(A[0, 0])
    if dtype != np.complex128:
        # the slice views work only properly for arrays
        A = np.asarray(A, dtype=float)

    # Quick return if possible
    if n % 2 == 1:
        return 0

    if not overwrite_a:
        A = A.copy()

    pfaffian_val = 1.0

    for k in range(0, n - 1, 2):
        # First, find the largest entry in A[k+1:,k] and
        # permute it to A[k+1,k]
        kp = k + 1 + np.abs(A[k + 1 :, k]).argmax()

        # Check if we need to pivot
        if kp != k + 1:
            # interchange rows k+1 and kp
            temp = A[k + 1, k:].copy()
            A[k + 1, k:] = A[kp, k:]
            A[kp, k:] = temp

            # Then interchange columns k+1 and kp
            temp = A[k:, k + 1].copy()
            A[k:, k + 1] = A[k:, kp]
            A[k:, kp] = temp

            # every interchange corresponds to a "-" in det(P)
            pfaffian_val *= -1

        # Now form the Gauss vector
        if A[k + 1, k] != 0.0:
            tau = A[k, k + 2 :].copy()
            tau = tau / A[k, k + 1]

            pfaffian_val *= A[k, k + 1]

            if k + 2 < n:
                # Update the matrix block A(k+2:,k+2)
                A[k + 2 :, k + 2 :] = A[k + 2 :, k + 2 :] + np.outer(
                    tau, A[k + 2 :, k + 1]
                )
                A[k + 2 :, k + 2 :] = A[k + 2 :, k + 2 :] - np.outer(
                    A[k + 2 :, k + 1], tau
                )
        else:
            # if we encounter a zero on the super/subdiagonal, the
            # Pfaffian is 0
            return 0.0

    return pfaffian_val

@cython.boundscheck(False)
@cython.wraparound(False)
def pfaffian_householder_real(cnp.ndarray[cnp.float64_t, ndim=2] A, bint overwrite_a=False):
    assert A.shape[0] == A.shape[1] and A.shape[0] > 0, "Matrix must be square and non-empty"
    assert fabs((A + A.T).max()) < 1e-14, "Matrix must be skew-symmetric"

    cdef:
        int n = A.shape[0], i
        cnp.float64_t pfaffian_val = 1.0
        cnp.float64_t alpha, tau
        cnp.ndarray[cnp.float64_t, ndim=1] v, w

    if n % 2 == 1:
        return 0.0

    if not overwrite_a:
        A = np.asarray(A.copy(), dtype=np.double)

    for i in range(n - 2):
        v = get_column_slice_as_1d_real(A, i + 1, i)
        v, tau, alpha = householder_real(v)
        A[i + 1, i] = alpha
        A[i, i + 1] = -alpha
        A[i + 2:, i] = 0
        A[i, i + 2:] = 0

        w = tau * np.dot(A[i + 1:, i + 1:], v)
        A[i + 1:, i + 1:] += np.outer(v, w) - np.outer(w, v)

        if tau != 0:
            pfaffian_val *= 1 - tau
        if i % 2 == 0:
            pfaffian_val *= -alpha

    pfaffian_val *= A[n - 2, n - 1]

    return pfaffian_val

@cython.boundscheck(False)
@cython.wraparound(False)
def pfaffian_householder_complex(cnp.ndarray[cnp.complex128_t, ndim=2] A, bint overwrite_a=False):
    assert A.shape[0] == A.shape[1] and A.shape[0] > 0, "Matrix must be square and non-empty"
    assert fabs((A + A.T).max()) < 1e-14, "Matrix must be skew-symmetric"

    cdef:
        int n = A.shape[0], i
        double complex pfaffian_val = 1.0 + 0j
        double complex alpha, tau
        cnp.ndarray[cnp.complex128_t, ndim=1] v, w

    if n % 2 == 1:
        return 0.0

    if not overwrite_a:
        A = np.asarray(A.copy(), dtype=np.complex128)

    for i in range(n - 2):
        v = get_column_slice_as_1d_complex(A, i + 1, i)
        v, tau, alpha = householder_complex(v)
        A[i + 1, i] = alpha
        A[i, i + 1] = -alpha
        A[i + 2:, i] = 0
        A[i, i + 2:] = 0

        w = tau * np.dot(A[i+1:, i+1:], v.conj())
        A[i + 1:, i + 1:] += np.outer(v, w) - np.outer(w, v)

        if tau != 0:
            pfaffian_val *= 1 - tau
        if i % 2 == 0:
            pfaffian_val *= -alpha

    pfaffian_val *= A[n - 2, n - 1]

    return pfaffian_val


@cython.boundscheck(False)
@cython.wraparound(False)
def pfaffian_householder(cnp.ndarray A, bint overwrite_a=False):
    # Type checks and setup
    assert A.shape[0] == A.shape[1] and A.shape[0] > 0, "Matrix must be square and non-empty"
    assert fabs((A + A.T).max()) < 1e-14, "Matrix must be skew-symmetric"

    # Dispatch to the appropriate function based on datatype
    if np.issubdtype(A.dtype, cnp.complexfloating):
        return pfaffian_householder_complex(A, overwrite_a=overwrite_a)
    elif np.issubdtype(A.dtype, np.number):
        return pfaffian_householder_real(A, overwrite_a=overwrite_a)
    else:
        raise TypeError("pfaffian_householder can only work on numeric input")

def pfaffian_schur(A, overwrite_a=False):
    """Calculate Pfaffian of a real antisymmetric matrix using
    the Schur decomposition. (Hessenberg would in principle be faster,
    but scipy-0.8 messed up the performance for scipy.linalg.hessenberg()).

    This function does not make use of the skew-symmetry of the matrix A,
    but uses a LAPACK routine that is coded in FORTRAN and hence faster
    than python. As a consequence, pfaffian_schur is only slightly slower
    than pfaffian().
    """

    assert np.issubdtype(A.dtype, np.number) and not np.issubdtype(
        A.dtype, np.complexfloating
    )

    assert A.shape[0] == A.shape[1] > 0

    assert abs(A + A.T).max() < 1e-14

    # Quick return if possible
    if A.shape[0] % 2 == 1:
        return 0

    (t, z) = la.schur(A, output="real", overwrite_a=overwrite_a)
    l = np.diag(t, 1)  # noqa: E741
    return np.prod(l[::2]) * la.det(z)
