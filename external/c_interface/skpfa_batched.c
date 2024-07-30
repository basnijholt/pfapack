#include "commondefs.h"
#include <stdio.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

int skpfa_batched_d(int batch_size, int N, double *A_batch, double *PFAFF_batch,
                    const char *UPLO, const char *MTHD)
{
    char uplo, mthd;

    uplo = toupper(UPLO[0]);
    mthd = toupper(MTHD[0]);

    if (N < 0) {
        return -1;
    }
    else if (A_batch == NULL) {
        return -2;
    }
    else if (PFAFF_batch == NULL) {
        return -3;
    }
    else if (uplo != 'U' && uplo != 'L') {
        return -4;
    }
    else if (mthd != 'P' && mthd != 'H') {
        return -5;
    }

    if (N > 0) {
        int ldim = N;
        int info = 0;
        int *iwork;
        double *work;
        int lwork = -1;
        double qwork;

        iwork = (int *)malloc(sizeof(int) * N);
        if (!iwork) return -100;

        // Workspace query
        PFAPACK_dskpfa(&uplo, &mthd, &N, A_batch, &ldim, &qwork,
                       iwork, &qwork, &lwork, &info);

        if (info) {
            free(iwork);
            return -101;
        }

        lwork = (int)qwork;
        work = (double *)malloc(sizeof(double) * lwork);
        if (!work) {
            free(iwork);
            return -102;
        }

        // Process each matrix in the batch
        for (int i = 0; i < batch_size; i++) {
            double *current_matrix = A_batch + i * N * N;
            double *current_pfaffian = PFAFF_batch + i;

            // Create a copy of the current matrix with the correct sign convention
            double *matrix_copy = (double *)malloc(sizeof(double) * N * N);
            if (!matrix_copy) {
                free(work);
                free(iwork);
                return -103;
            }

            // Copy the upper triangular part and negate it
            for (int row = 0; row < N; row++) {
                for (int col = row + 1; col < N; col++) {
                    matrix_copy[row * N + col] = -current_matrix[row * N + col];
                    matrix_copy[col * N + row] = current_matrix[row * N + col];
                }
                matrix_copy[row * N + row] = 0.0;  // Diagonal elements should be zero
            }

            PFAPACK_dskpfa(&uplo, &mthd, &N, matrix_copy, &ldim, current_pfaffian,
                           iwork, work, &lwork, &info);

            free(matrix_copy);

            if (info) {
                free(work);
                free(iwork);
                return -104 - i;
            }
        }

        free(work);
        free(iwork);
    }
    else {
        for (int i = 0; i < batch_size; i++) {
            PFAFF_batch[i] = 1.0;
        }
    }

    return 0;
}

int skpfa_batched_z(int batch_size, int N, double *A_batch_real_imag, double *PFAFF_batch_real_imag,
                    const char *UPLO, const char *MTHD)
{
    char uplo, mthd;

    uplo = toupper(UPLO[0]);
    mthd = toupper(MTHD[0]);

    if (N < 0) {
        return -1;
    }
    else if (A_batch_real_imag == NULL) {
        return -2;
    }
    else if (PFAFF_batch_real_imag == NULL) {
        return -3;
    }
    else if (uplo != 'U' && uplo != 'L') {
        return -4;
    }
    else if (mthd != 'P' && mthd != 'H') {
        return -5;
    }

    if (N > 0) {
        int ldim = N;
        int info = 0;
        int *iwork, lwork;
        double *rwork;
        double complex *work, qwork;

        iwork = (int *)malloc(sizeof(int) * N);
        if (!iwork) {
            return -100;
        }
        rwork = (double *)malloc(sizeof(double) * (N-1));
        if (!rwork) {
            free(iwork);
            return -100;
        }

        // Workspace query
        lwork = -1;
        double complex dummy_matrix[1];
        double complex dummy_pfaffian;
        PFAPACK_zskpfa(&uplo, &mthd, &N, dummy_matrix, &ldim, &dummy_pfaffian,
                       iwork, &qwork, &lwork, rwork, &info);

        if (info) {
            free(rwork);
            free(iwork);
            return -101;
        }

        lwork = (int)creal(qwork);

        work = (double complex *)malloc(sizeof(double complex) * lwork);
        if (!work) {
            // Try minimal workspace
            if (mthd == 'P') lwork = 1;
            else lwork = 2*N-1;

            work = (double complex *)malloc(sizeof(double complex) * lwork);
            if (!work) {
                free(rwork);
                free(iwork);
                return -100;
            }
        }

        // Process each matrix in the batch
        for (int i = 0; i < batch_size; i++) {
            double *current_matrix_real = A_batch_real_imag + i * 2 * N * N;
            double *current_matrix_imag = current_matrix_real + N * N;
            double *current_pfaffian_real = PFAFF_batch_real_imag + i * 2;
            double *current_pfaffian_imag = current_pfaffian_real + 1;

            double complex *matrix_copy = (double complex *)malloc(sizeof(double complex) * N * N);
              if (!matrix_copy) {
                  free(work);
                  free(rwork);
                  free(iwork);
                  return -103;
              }
              
              // Copy the upper triangular part and negate it
              for (int row = 0; row < N; row++) {
                  for (int col = row + 1; col < N; col++) {
                      double complex value = current_matrix_real[row * N + col] + I * current_matrix_imag[row * N + col];
                      matrix_copy[row * N + col] = -value;
                      matrix_copy[col * N + row] = value;
                  }
                  matrix_copy[row * N + row] = 0.0;  // Diagonal elements should be zero
              }

            double complex current_pfaffian;
            PFAPACK_zskpfa(&uplo, &mthd, &N, matrix_copy, &ldim, &current_pfaffian,
                           iwork, work, &lwork, rwork, &info);

            *current_pfaffian_real = creal(current_pfaffian);
            *current_pfaffian_imag = cimag(current_pfaffian);

            free(matrix_copy);

            if (info) {
                free(work);
                free(rwork);
                free(iwork);
                return -104 - i;
            }
        }

        free(work);
        free(rwork);
        free(iwork);
    }
    else {
        for (int i = 0; i < batch_size; i++) {
            PFAFF_batch_real_imag[2*i] = 1.0;
            PFAFF_batch_real_imag[2*i+1] = 0.0;
        }
    }

    return 0;
}

#ifdef __cplusplus
}
#endif
