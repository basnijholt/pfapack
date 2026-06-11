/**********************************************************************
  Purpose
  =======

  skpfa_batched_x computes the Pfaffians of a batch of dense
  skew-symmetric matrices. Compared to calling skpfa_x in a loop, the
  LAPACK workspace is allocated once for the whole batch rather than
  once per matrix, which significantly reduces overhead for batches of
  many small matrices.

  Usage
  =====

  int skpfa_batched_x(int BATCH_SIZE, int N, const datatype *A_BATCH,
                      datatype *PFAFF_BATCH, const char *UPLO,
                      const char *MTHD)

  where
     datatype = double or doublecmplx
  for x = d or z

  Arguments
  =========

  BATCH_SIZE (input) int
          Number of matrices in the batch. BATCH_SIZE >= 0.

  N       (input) int
          Size of each matrix. N >= 0.

  A_BATCH (input) datatype *
          pointer to a memory block of size
          BATCH_SIZE*N*N*sizeof(datatype) holding BATCH_SIZE
          consecutive skew-symmetric NxN-matrices in C format
          (row-major order). Note that this differs from skpfa_x,
          which expects Fortran format. The input is not modified.
             If UPLO = 'U', the upper triangular part of each matrix
                contains its upper triangular part, and the strictly
                lower triangular part is not referenced.
             If UPLO = 'L', vice versa.

  PFAFF_BATCH (output) datatype *
          pointer to a memory block of size BATCH_SIZE*sizeof(datatype)
          The values of the Pfaffians.

  UPLO    (input) char *
            = 'U':  Upper triangle of each matrix is stored;
            = 'L':  Lower triangle of each matrix is stored.

  MTHD    (input) char *
            = 'P': Compute Pfaffian using Parlett-Reid algorithm
                   (recommended)
            = 'H': Compute Pfaffian using Householder reflections

  Return value
  ============

  The return value indicates whether an error occured:
      0:    successful exit
    < 0:    if the return value is -i, the i-th argument had an
            illegal value
   -100:    failed to allocate enough internal memory

  Further Details
  ===============

  The row-major matrices are fed directly to the column-major Fortran
  kernels without transposing: the kernel then operates on A^T rather
  than A. This is accounted for by swapping UPLO and using the identity
  pf(A^T) = pf(-A) = (-1)^(N/2) pf(A) to correct the sign of the
  result, avoiding an explicit transpose of every matrix in the batch.

  For small matrices (N <= PF_SMALL_N_MAX) and MTHD = 'P', a hand-rolled
  Parlett-Reid elimination with partial pivoting is used instead of the
  Fortran kernels. It performs the same algorithm with the same pivoting
  strategy, but avoids the LAPACK calling machinery (argument checking,
  workspace logic, BLAS calls), which dominates the runtime at small N.

  Ported from the batched implementation by Joshua Goings (@jjgoings),
  https://github.com/jjgoings/pfapack.

****************************************************************************/

#include <math.h>
#include <string.h>

#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Below this size the hand-rolled Parlett-Reid kernels beat the Fortran
   path; above it LAPACK's blocked algorithms win. */
#define PF_SMALL_N_MAX 24

#define A(i, j) a[(size_t)(i) * n + (j)]

/* Expand the UPLO triangle of a row-major skew-symmetric matrix into a
   full scratch matrix, referencing only that triangle. Writes are kept
   linear; the mirrored values are read back from the scratch rows
   already filled, which stay in L1 at these sizes. */
static void pf_small_fill_d(double *a, const double *src, int n, char uplo)
{
    int i, j;
    if (uplo == 'U') {
        for (i = 0; i < n; i++) {
            const double *srow = src + (size_t)i * n;
            for (j = 0; j < i; j++) A(i, j) = -A(j, i);
            A(i, i) = 0.0;
            for (j = i + 1; j < n; j++) A(i, j) = srow[j];
        }
    } else {
        for (i = 0; i < n; i++) {
            const double *srow = src + (size_t)i * n;
            for (j = 0; j < i; j++) A(i, j) = srow[j];
            A(i, i) = 0.0;
            for (j = i + 1; j < n; j++) A(i, j) = -src[(size_t)j * n + i];
        }
    }
}

/* Parlett-Reid LTL^T elimination with partial pivoting on a full
   row-major skew-symmetric matrix (destroyed). Mirrors pfaffian_LTL
   from pfapack/pfaffian.py. N must be even. */
static double pf_small_d(double *a, int n)
{
    double tau[PF_SMALL_N_MAX], col[PF_SMALL_N_MAX];
    double pfaff = 1.0;
    int i, j, k;

    for (k = 0; k + 1 < n; k += 2) {
        /* pivot: largest |A[k+1:, k]| */
        int kp = k + 1;
        double mx = fabs(A(k + 1, k));
        for (i = k + 2; i < n; i++) {
            double v = fabs(A(i, k));
            if (v > mx) { mx = v; kp = i; }
        }
        if (kp != k + 1) {
            for (j = k; j < n; j++) {
                double t = A(k + 1, j); A(k + 1, j) = A(kp, j); A(kp, j) = t;
            }
            for (i = k; i < n; i++) {
                double t = A(i, k + 1); A(i, k + 1) = A(i, kp); A(i, kp) = t;
            }
            pfaff = -pfaff;
        }
        if (A(k + 1, k) == 0.0) return 0.0;
        pfaff *= A(k, k + 1);
        if (k + 2 < n) {
            /* rank-2 update of the trailing block:
               A[k+2:,k+2:] += tau (x) col - col (x) tau,
               tau = A[k, k+2:] / A[k, k+1], col = A[k+2:, k+1] */
            double inv_piv = 1.0 / A(k, k + 1);
            for (j = k + 2; j < n; j++) {
                tau[j] = A(k, j) * inv_piv;
                col[j] = A(j, k + 1);
            }
            for (i = k + 2; i < n; i++) {
                double ti = tau[i], ci = col[i];
                double *row = &A(i, 0);
                for (j = k + 2; j < n; j++)
                    row[j] += ti * col[j] - ci * tau[j];
            }
        }
    }
    return pfaff;
}

/* Complex twin of the above. Needs native complex arithmetic, so it is
   only available with C99 or C++ complex; the STRUCT_COMPLEX fallback
   keeps using the Fortran path. */
#if defined(C99_COMPLEX) || defined(CPLUSPLUS_COMPLEX)
#define PF_HAVE_SMALL_Z 1

static double pf_small_cabs1(doublecmplx x)
{
    return fabs(real(x)) + fabs(imag(x));
}

static void pf_small_fill_z(doublecmplx *a, const doublecmplx *src, int n,
                            char uplo)
{
    int i, j;
    if (uplo == 'U') {
        for (i = 0; i < n; i++) {
            const doublecmplx *srow = src + (size_t)i * n;
            for (j = 0; j < i; j++) A(i, j) = -A(j, i);
            A(i, i) = doublecmplx_zero;
            for (j = i + 1; j < n; j++) A(i, j) = srow[j];
        }
    } else {
        for (i = 0; i < n; i++) {
            const doublecmplx *srow = src + (size_t)i * n;
            for (j = 0; j < i; j++) A(i, j) = srow[j];
            A(i, i) = doublecmplx_zero;
            for (j = i + 1; j < n; j++) A(i, j) = -src[(size_t)j * n + i];
        }
    }
}

static doublecmplx pf_small_z(doublecmplx *a, int n)
{
    doublecmplx tau[PF_SMALL_N_MAX], col[PF_SMALL_N_MAX];
    doublecmplx pfaff = doublecmplx_one;
    int i, j, k;

    for (k = 0; k + 1 < n; k += 2) {
        int kp = k + 1;
        double mx = pf_small_cabs1(A(k + 1, k));
        for (i = k + 2; i < n; i++) {
            double v = pf_small_cabs1(A(i, k));
            if (v > mx) { mx = v; kp = i; }
        }
        if (kp != k + 1) {
            for (j = k; j < n; j++) {
                doublecmplx t = A(k + 1, j);
                A(k + 1, j) = A(kp, j);
                A(kp, j) = t;
            }
            for (i = k; i < n; i++) {
                doublecmplx t = A(i, k + 1);
                A(i, k + 1) = A(i, kp);
                A(i, kp) = t;
            }
            pfaff = -pfaff;
        }
        if (mx == 0.0) return doublecmplx_zero;
        pfaff *= A(k, k + 1);
        if (k + 2 < n) {
            doublecmplx inv_piv = doublecmplx_one / A(k, k + 1);
            for (j = k + 2; j < n; j++) {
                tau[j] = A(k, j) * inv_piv;
                col[j] = A(j, k + 1);
            }
            for (i = k + 2; i < n; i++) {
                doublecmplx ti = tau[i], ci = col[i];
                doublecmplx *row = &A(i, 0);
                for (j = k + 2; j < n; j++)
                    row[j] += ti * col[j] - ci * tau[j];
            }
        }
    }
    return pfaff;
}
#endif /* C99_COMPLEX || CPLUSPLUS_COMPLEX */

#undef A

int skpfa_batched_d(int batch_size, int N,
                    const double *A_batch, double *PFAFF_batch,
                    const char *UPLO, const char *MTHD)
{
    char uplo = toupper(UPLO[0]);
    char mthd = toupper(MTHD[0]);
    int i, info = 0, ret = 0;

    if (batch_size < 0) return -1;
    if (N < 0) return -2;
    if (A_batch == NULL) return -3;
    if (PFAFF_batch == NULL) return -4;
    if (uplo != 'U' && uplo != 'L') return -5;
    if (mthd != 'P' && mthd != 'H') return -6;

    if (N == 0) {
        for (i = 0; i < batch_size; i++) PFAFF_batch[i] = 1.0;
        return 0;
    }

    if (N % 2 != 0) {
        for (i = 0; i < batch_size; i++) PFAFF_batch[i] = 0.0;
        return 0;
    }

    if (N <= PF_SMALL_N_MAX && mthd == 'P') {
        double scratch[PF_SMALL_N_MAX * PF_SMALL_N_MAX];
        const size_t matrix_elems = (size_t)N * (size_t)N;
        for (i = 0; i < batch_size; i++) {
            pf_small_fill_d(scratch, A_batch + (size_t)i * matrix_elems,
                            N, uplo);
            PFAFF_batch[i] = pf_small_d(scratch, N);
        }
        return 0;
    }

    /* Row-major input read as column-major is A^T, whose data lives in
       the opposite triangle; the sign is fixed up after the fact. */
    {
        const char uplo_t = (uplo == 'U') ? 'L' : 'U';
        const int sign_flip = (N / 2) % 2;
        const size_t matrix_elems = (size_t)N * (size_t)N;
        int ldim = N, lwork = -1;
        double qwork = 0.0, pfaff_dummy = 0.0;
        double *work;

        double *A_work = (double *)malloc(matrix_elems * sizeof(double));
        int *iwork = (int *)malloc((size_t)N * sizeof(int));
        if (!A_work || !iwork) {
            free(iwork);
            free(A_work);
            return -100;
        }

        /* Workspace query, done once for the whole batch */
        PFAPACK_dskpfa(&uplo_t, &mthd, &N, A_work, &ldim, &pfaff_dummy,
                       iwork, &qwork, &lwork, &info);
        if (info) {
            free(iwork);
            free(A_work);
            return info;
        }

        lwork = (int)qwork;
        if (lwork < 1) lwork = 1;
        work = (double *)malloc((size_t)lwork * sizeof(double));
        if (!work) {
            /* Fall back to the minimal workspace */
            lwork = (mthd == 'P') ? 1 : (2 * N - 1);
            work = (double *)malloc((size_t)lwork * sizeof(double));
            if (!work) {
                free(iwork);
                free(A_work);
                return -100;
            }
        }

        for (i = 0; i < batch_size; i++) {
            memcpy(A_work, A_batch + (size_t)i * matrix_elems,
                   matrix_elems * sizeof(double));
            PFAPACK_dskpfa(&uplo_t, &mthd, &N, A_work, &ldim,
                           &PFAFF_batch[i], iwork, work, &lwork, &info);
            if (info) {
                ret = info;
                break;
            }
            if (sign_flip) PFAFF_batch[i] = -PFAFF_batch[i];
        }

        free(work);
        free(iwork);
        free(A_work);
    }

    return ret;
}

int skpfa_batched_z(int batch_size, int N,
                    const doublecmplx *A_batch, doublecmplx *PFAFF_batch,
                    const char *UPLO, const char *MTHD)
{
    char uplo = toupper(UPLO[0]);
    char mthd = toupper(MTHD[0]);
    int i, info = 0, ret = 0;

    if (batch_size < 0) return -1;
    if (N < 0) return -2;
    if (A_batch == NULL) return -3;
    if (PFAFF_batch == NULL) return -4;
    if (uplo != 'U' && uplo != 'L') return -5;
    if (mthd != 'P' && mthd != 'H') return -6;

    if (N == 0) {
        for (i = 0; i < batch_size; i++) PFAFF_batch[i] = doublecmplx_one;
        return 0;
    }

    if (N % 2 != 0) {
        for (i = 0; i < batch_size; i++) PFAFF_batch[i] = doublecmplx_zero;
        return 0;
    }

#ifdef PF_HAVE_SMALL_Z
    if (N <= PF_SMALL_N_MAX && mthd == 'P') {
        doublecmplx scratch[PF_SMALL_N_MAX * PF_SMALL_N_MAX];
        const size_t matrix_elems = (size_t)N * (size_t)N;
        for (i = 0; i < batch_size; i++) {
            pf_small_fill_z(scratch, A_batch + (size_t)i * matrix_elems,
                            N, uplo);
            PFAFF_batch[i] = pf_small_z(scratch, N);
        }
        return 0;
    }
#endif

    /* Row-major input read as column-major is A^T, whose data lives in
       the opposite triangle; the sign is fixed up after the fact. */
    {
        const char uplo_t = (uplo == 'U') ? 'L' : 'U';
        const int sign_flip = (N / 2) % 2;
        const size_t matrix_elems = (size_t)N * (size_t)N;
        int ldim = N, lwork = -1;
        doublecmplx qwork = doublecmplx_zero;
        doublecmplx pfaff_dummy = doublecmplx_zero;
        doublecmplx *work;
        double *rwork;

        doublecmplx *A_work =
            (doublecmplx *)malloc(matrix_elems * sizeof(doublecmplx));
        int *iwork = (int *)malloc((size_t)N * sizeof(int));
        rwork = (double *)malloc((size_t)(N > 1 ? N - 1 : 1) * sizeof(double));
        if (!A_work || !iwork || !rwork) {
            free(rwork);
            free(iwork);
            free(A_work);
            return -100;
        }

        /* Workspace query, done once for the whole batch */
        PFAPACK_zskpfa(&uplo_t, &mthd, &N, A_work, &ldim, &pfaff_dummy,
                       iwork, &qwork, &lwork, rwork, &info);
        if (info) {
            free(rwork);
            free(iwork);
            free(A_work);
            return info;
        }

        lwork = (int)real(qwork);
        if (lwork < 1) lwork = 1;
        work = (doublecmplx *)malloc((size_t)lwork * sizeof(doublecmplx));
        if (!work) {
            /* Fall back to the minimal workspace */
            lwork = (mthd == 'P') ? 1 : (2 * N - 1);
            work = (doublecmplx *)malloc((size_t)lwork * sizeof(doublecmplx));
            if (!work) {
                free(rwork);
                free(iwork);
                free(A_work);
                return -100;
            }
        }

        for (i = 0; i < batch_size; i++) {
            memcpy(A_work, A_batch + (size_t)i * matrix_elems,
                   matrix_elems * sizeof(doublecmplx));
            PFAPACK_zskpfa(&uplo_t, &mthd, &N, A_work, &ldim,
                           &PFAFF_batch[i], iwork, work, &lwork, rwork,
                           &info);
            if (info) {
                ret = info;
                break;
            }
#ifdef STRUCT_COMPLEX
            if (sign_flip) {
                PFAFF_batch[i].re = -PFAFF_batch[i].re;
                PFAFF_batch[i].im = -PFAFF_batch[i].im;
            }
#else
            if (sign_flip) PFAFF_batch[i] = -PFAFF_batch[i];
#endif
        }

        free(work);
        free(rwork);
        free(iwork);
        free(A_work);
    }

    return ret;
}

#ifdef __cplusplus
}
#endif
