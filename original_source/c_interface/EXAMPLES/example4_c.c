#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pfapack.h"

/* matrix multiplication routine from BLAS */
#define BLAS_dgemm fortran_name(dgemm, DGEMM)
/* assemble Householder reflections, from LAPACK */

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
  /* dense real example */
  int N=4;
  int info, i, j;
  double one=1, zero=0;

  double *A=(double *)malloc(sizeof(double)*N*N);
  double *L=(double *)malloc(sizeof(double)*N*N);
  double *T=(double *)malloc(sizeof(double)*N*N);
  int *ipiv=(int *)malloc(sizeof(int)*N);

  /* build up a skewsymmetric matrix */
  memset(A, 0, sizeof(double)*N*N);

  A[dense_fortran(1,2,N)]=1.0;
  A[dense_fortran(2,1,N)]=-1.0;

  A[dense_fortran(1,3,N)]=2.0;
  A[dense_fortran(3,1,N)]=-2.0;

  A[dense_fortran(1,4,N)]=3.0;
  A[dense_fortran(4,1,N)]=-3.0;

  A[dense_fortran(2,3,N)]=4.0;
  A[dense_fortran(3,2,N)]=-4.0;

  A[dense_fortran(2,4,N)]=5.0;
  A[dense_fortran(4,2,N)]=-5.0;

  A[dense_fortran(3,4,N)]=6.0;
  A[dense_fortran(4,3,N)]=-6.0;

  printf("The original matrix was:\n\n");
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 A[dense_fortran(1,1,N)], A[dense_fortran(1,2,N)],
	 A[dense_fortran(1,3,N)], A[dense_fortran(1,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 A[dense_fortran(2,1,N)], A[dense_fortran(2,2,N)],
	 A[dense_fortran(2,3,N)], A[dense_fortran(2,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 A[dense_fortran(3,1,N)], A[dense_fortran(3,2,N)],
	 A[dense_fortran(3,3,N)], A[dense_fortran(3,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 A[dense_fortran(4,1,N)], A[dense_fortran(4,2,N)],
	 A[dense_fortran(4,3,N)], A[dense_fortran(4,4,N)]);

  /* compute the L T L^T decomposition */

  info = sktrf_d(N, A, ipiv, "L", "N");
  assert(info == 0);

  printf("\nThe matrix P A P^T is factorized as L * T * L^T with\n");
  printf("tridiagonal matrix T:\n\n");
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 0.0, A[dense_fortran(1,2,N)], 0.0, 0.0);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 -A[dense_fortran(1,2,N)], 0.0, A[dense_fortran(2,3,N)], 0.0);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 0.0, -A[dense_fortran(2,3,N)], 0.0, A[dense_fortran(3,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 0.0, 0.0, -A[dense_fortran(3,4,N)], 0.0);

  /* Copy the tridiagonal matrix T to dense format for easier handling */
  memset(T, 0, sizeof(double) *N*N);

  for(i=1; i<N; i++) {
    T[dense_fortran(i+1,i,N)] = A[dense_fortran(i+1,i,N)];
    T[dense_fortran(i,i+1,N)] = -T[dense_fortran(i+1,i,N)];
  }

  /* Copy the lower triangular matrix L to dense format for easier handling */
  memset(L, 0, sizeof(double) *N*N);

  for(j=1; j<=N; j++) {
    L[dense_fortran(j,j,N)]=1;

    for(i=j+2; i<=N; i++) {
      L[dense_fortran(i,j+1,N)] = A[dense_fortran(i,j,N)];
    }
  }

  printf("\nlower triangular matrix L:\n\n");
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 L[dense_fortran(1,1,N)], L[dense_fortran(1,2,N)],
	 L[dense_fortran(1,3,N)], L[dense_fortran(1,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 L[dense_fortran(2,1,N)], L[dense_fortran(2,2,N)],
	 L[dense_fortran(2,3,N)], L[dense_fortran(2,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 L[dense_fortran(3,1,N)], L[dense_fortran(3,2,N)],
	 L[dense_fortran(3,3,N)], L[dense_fortran(3,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 L[dense_fortran(4,1,N)], L[dense_fortran(4,2,N)],
	 L[dense_fortran(4,3,N)], L[dense_fortran(4,4,N)]);


  printf("\nSanity check: P^T * L * T * L^T *P should give (approximately) the original matrix:\n\n");

  /* A = (L * T) */
  BLAS_dgemm("N", "N", &N, &N, &N, &one, L, &N, T, &N, &zero, A, &N);

  /* T = (A * L^T) */
  BLAS_dgemm("N", "T", &N, &N, &N, &one, A, &N, L, &N, &zero, T, &N);

  /* apply the inverse permutation to L*T*L^T
     This is done by doing the interchanges described in IPIV from N=4 to 1
     (we need to apply the *inverse* permutation)
  */
  for(i=N; i>=1; i--) {
    /* interchange rows first */
    for(j=1; j<=N; j++) {
      double tmp;

      tmp = T[dense_fortran(i,j,N)];
      T[dense_fortran(i,j,N)] = T[dense_fortran(ipiv[i-1],j,N)];
      T[dense_fortran(ipiv[i-1],j,N)] = tmp;
    }

    /* then interchange columns */
    for(j=1; j<=N; j++) {
      double tmp;

      tmp = T[dense_fortran(j,i,N)];
      T[dense_fortran(j,i,N)] = T[dense_fortran(j,ipiv[i-1],N)];
      T[dense_fortran(j,ipiv[i-1],N)] = tmp;
    }
  }

  printf("%8.2g %8.2g %8.2g %8.2g\n",
	 T[dense_fortran(1,1,N)], T[dense_fortran(1,2,N)],
	 T[dense_fortran(1,3,N)], T[dense_fortran(1,4,N)]);
  printf("%8.2g %8.2g %8.2g %8.2g\n",
	 T[dense_fortran(2,1,N)], T[dense_fortran(2,2,N)],
	 T[dense_fortran(2,3,N)], T[dense_fortran(2,4,N)]);
  printf("%8.2g %8.2g %8.2g %8.2g\n",
	 T[dense_fortran(3,1,N)], T[dense_fortran(3,2,N)],
	 T[dense_fortran(3,3,N)], T[dense_fortran(3,4,N)]);
  printf("%8.2g %8.2g %8.2g %8.2g\n",
	 T[dense_fortran(4,1,N)], T[dense_fortran(4,2,N)],
	 T[dense_fortran(4,3,N)], T[dense_fortran(4,4,N)]);

  return 0;
}
