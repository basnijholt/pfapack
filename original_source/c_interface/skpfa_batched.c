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

  Ported from the batched implementation by Joshua Goings (@jjgoings),
  https://github.com/jjgoings/pfapack.

****************************************************************************/

#include <string.h>

#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

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
