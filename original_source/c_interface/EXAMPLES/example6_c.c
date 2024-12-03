#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pfapack.h"

/* matrix multiplication routine from BLAS */
#define BLAS_dgemm fortran_name(dgemm, DGEMM)
#ifdef __cplusplus
extern "C" {
#endif
void BLAS_dgemm(const char *, const char *, const int *,
		const int *, const int *, const double *,
		const double *, const int *, const double *,
		const int *, const double *, double *,
		const int *);
#ifdef __cplusplus
}
#endif

int main()
{
  /* banded real example */
  int N=4;
  int KD=2;
  int KD1=KD+1;
  int info, i;
  double detQ;
  double one = 1, zero = 0;

  double *Ad = (double *)malloc(sizeof(double)* KD1*N);
  double *Q = (double *)malloc(sizeof(double)* N*N);
  double *T = (double *)malloc(sizeof(double)* N*N);
  double *tmp = (double *)malloc(sizeof(double)* N*N);

  /* build up a skewsymmetric band matrix */

  memset(Ad, 0, sizeof(double)* KD1*N);

  Ad[bandlower_fortran(3,1,N,KD)]=1.0;

  Ad[bandlower_fortran(4,2,N,KD)]=1.0;

  /* compute the tridiagonal form under orthogonal similarity */

  info = skbtrd_d(N, KD, Ad, &detQ, Q, "L", "N", "V");
  assert(info == 0);

  printf("The original matrix was:\n\n 0  0 -1  0\n 0  0  0 -1\n 1  0  0  0\n 0  1  0  0\n\n");

  printf("\nThe matrix A is factorized as Q * T * Q^T with\n");
  printf("tridiagonal matrix T:\n\n");
  printf("%2g %2g %2g %2g\n",
	 0.0, -Ad[bandlower_fortran(2,1,N,KD)], 0.0, 0.0);
  printf("%2g %2g %2g %2g\n",
	 Ad[bandlower_fortran(2,1,N,KD)], 0.0, -Ad[bandlower_fortran(3,2,N,KD)], 0.0);
  printf("%2g %2g %2g %2g\n",
	 0.0, Ad[bandlower_fortran(3,2,N,KD)], 0.0, -Ad[bandlower_fortran(4,3,N,KD)]);
  printf("%2g %2g %2g %2g\n",
	 0.0, 0.0, Ad[bandlower_fortran(4,3,N,KD)], 0.0);

  printf("\nand orthogonal similarity transform Q:\n\n");
  printf("%2g %2g %2g %2g\n",
	 Q[dense_fortran(1,1,N)], Q[dense_fortran(1,2,N)],
	 Q[dense_fortran(1,3,N)], Q[dense_fortran(1,4,N)]);
  printf("%2g %2g %2g %2g\n",
	 Q[dense_fortran(2,1,N)], Q[dense_fortran(2,2,N)],
	 Q[dense_fortran(2,3,N)], Q[dense_fortran(2,4,N)]);
  printf("%2g %2g %2g %2g\n",
	 Q[dense_fortran(3,1,N)], Q[dense_fortran(3,2,N)],
	 Q[dense_fortran(3,3,N)], Q[dense_fortran(3,4,N)]);
  printf("%2g %2g %2g %2g\n",
	 Q[dense_fortran(4,1,N)], Q[dense_fortran(4,2,N)],
	 Q[dense_fortran(4,3,N)], Q[dense_fortran(4,4,N)]);

  printf("\nSanity check: Q * T * Q^T should give the original matrix:\n\n");

  /* Copy the tridiagonal matrix to dense format for easier handling */
  memset(T, 0, sizeof(double) *N*N);

  for(i=1; i<N; i++) {
    T[dense_fortran(i+1,i,N)] = Ad[bandlower_fortran(i+1,i,N,KD)];
    T[dense_fortran(i,i+1,N)] = -T[dense_fortran(i+1,i,N)];
  }

  /* tmp = (Q * T) */
  BLAS_dgemm("N", "N", &N, &N, &N, &one, Q, &N, T, &N, &zero, tmp, &N);

  /* T = (tmp * Q^T) */
  BLAS_dgemm("N", "T", &N, &N, &N, &one, tmp, &N, Q, &N, &zero, T, &N);

  printf("%2g %2g %2g %2g\n",
	 T[dense_fortran(1,1,N)], T[dense_fortran(1,2,N)],
	 T[dense_fortran(1,3,N)], T[dense_fortran(1,4,N)]);
  printf("%2g %2g %2g %2g\n",
	 T[dense_fortran(2,1,N)], T[dense_fortran(2,2,N)],
	 T[dense_fortran(2,3,N)], T[dense_fortran(2,4,N)]);
  printf("%2g %2g %2g %2g\n",
	 T[dense_fortran(3,1,N)], T[dense_fortran(3,2,N)],
	 T[dense_fortran(3,3,N)], T[dense_fortran(3,4,N)]);
  printf("%2g %2g %2g %2g\n",
	 T[dense_fortran(4,1,N)], T[dense_fortran(4,2,N)],
	 T[dense_fortran(4,3,N)], T[dense_fortran(4,4,N)]);

  return 0;
}
