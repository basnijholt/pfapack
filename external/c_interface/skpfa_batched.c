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

    // ... (error checking code omitted for brevity)

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

#ifdef __cplusplus
}
#endif
