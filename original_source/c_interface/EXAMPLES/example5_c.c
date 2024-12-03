#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pfapack.h"

/* matrix multiplication routine from BLAS */
#define BLAS_dgemm fortran_name(dgemm, DGEMM)
/* assemble Householder reflections, from LAPACK */
#define LAPACK_dorgtr fortran_name(dorgtr, DORGTR)

#ifdef __cplusplus
extern "C" {
#endif
void BLAS_dgemm(const char *, const char *, const int *,
		const int *, const int *, const double *,
		const double *, const int *, const double *,
		const int *, const double *, double *,
		const int *);

void LAPACK_dorgtr(const char *, const int *, double *,
		   const int *, const double *, double *,
		   const int *, int *);
#ifdef __cplusplus
}
#endif

int main()
{
  /* dense real example */
  int N=4;
  int info, i;
  int lwork=N-1;
  double one=1, zero=0;

  double *A=(double *)malloc(sizeof(double)*N*N);
  double *tmp=(double *)malloc(sizeof(double)*N*N);
  double *T=(double *)malloc(sizeof(double)*N*N);
  double *tau=(double *)malloc(sizeof(double)*(N-1));
  double *work=(double *)malloc(sizeof(double)*(N-1));

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

  /* compute the tridiagonal form under orthogonal similarity */

  info = sktrd_d(N, A, tau, "U", "N");
  assert(info == 0);

  printf("\nThe matrix A is factorized as Q * T * Q^T with\n");
  printf("tridiagonal matrix T:\n\n");
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 0.0, A[dense_fortran(1,2,N)], 0.0, 0.0);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 -A[dense_fortran(1,2,N)], 0.0, A[dense_fortran(2,3,N)], 0.0);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 0.0, -A[dense_fortran(2,3,N)], 0.0, A[dense_fortran(3,4,N)]);
  printf("%5.2g %5.2g %5.2g %5.2g\n",
	 0.0, 0.0, -A[dense_fortran(3,4,N)], 0.0);

  /* Copy the tridiagonal matrix to dense format for easier handling */
  /* (A is overwritten by Q in the next step) */
  memset(T, 0, sizeof(double) *N*N);

  for(i=1; i<N; i++) {
    T[dense_fortran(i,i+1,N)] = A[dense_fortran(i,i+1,N)];
    T[dense_fortran(i+1,i,N)] = -T[dense_fortran(i,i+1,N)];
  }

  /* form the unitary matrix Q from the Householder reflections */
  /* (call uses only minimal workspace) */
  LAPACK_dorgtr("U", &N, A, &N, tau, work, &lwork, &info);

  printf("\nand orthogonal similarity transform Q:\n\n");
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

  printf("\nSanity check: Q * T * Q^T should give (approximately) the original matrix:\n\n");

  /* tmp = (Q * T) */
  BLAS_dgemm("N", "N", &N, &N, &N, &one, A, &N, T, &N, &zero, tmp, &N);

  /* T = (tmp * Q^T) */
  BLAS_dgemm("N", "T", &N, &N, &N, &one, tmp, &N, A, &N, &zero, T, &N);

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
