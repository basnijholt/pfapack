/**********************************************************************
  Purpose
  =======

  sktrf_x computes the factorization of a skew-symmetric matrix A
  using the Parlett-Reid algorithm:

     P*A*P^T = U*T*U^T  or  P*A*P^T = L*T*L^T

  where U (or L) unit upper (lower) triangular matrix (^T denotes
  the transpose), T is a skew-symmetric tridiagonal matrix and P
  is a permutation matrix. In addition to being unit triangular,
  U(1:n-1,n)=0 and L(2:n,1)=0.
  Instead of a full tridiagonalization, SKTRF can also compute a
  partial tridiagonal form for computing the Pfaffian.

  This subroutine uses, if enough memory is available, the
  blocked version of the algorithm.

  Usage
  =====

  int sktrf_x(int N, datatype *A, int *IPIV,
	      const char *UPLO, const char *MODE)

  where
     datatype = float, double, floatcmplx, or doublecmplx
  for x = s, d, c, or z

  When using a C++ compiler, the trailing _x may be omitted.

  Arguments
  =========

  N       (input) int
          Size of the matrix A. If MODE = 'P', N must be even. N >= 0.

  A       (input/output) datatype *
          pointer to a memory block of size N*N*sizeof(datatype)
          On entry, the skew-symmetric NxN-matrix A in Fortran format.
             If UPLO = 'U', the upper triangular part of A contains
                the upper triangular part of the matrix A, and the
                strictly lower triangular part of A is not referenced.
             If UPLO = 'L', the lower triangular part of A contains
                the lower triangular part of the matrix A, and the
                strictly upper triangular part of A is not referenced.
          On exit, the tridiagonal matrix T and the multipliers used
          to obtain the factor U or L (see below for further details).

  IPIV    (output) int *,
          pointer to array of dimension N
          Information about the permutation matrix P: row and column
          i are interchanged with IPIV[i-1]. If UPLO = 'U', those
          interchanges are done in the order i = N ... 1, if UPLO = 'L'
          in the order i = 1 ... N.

  UPLO    (input) char *
            = 'U':  Upper triangular
            = 'L':  Lower triangular

  MODE    (input) char *
            = 'N':  A is fully tridiagonalized
            = 'P':  A is partially tridiagonalized for Pfaffian computation
                   (details see below)

  Return value
  ============

  The return value of sktrf_x indicates whether an error occured:
    >=0:    successful exit
            If the return value is k>0, then the off-diagonal entry in
            the k-th row (UPLO = 'U') or k-th column (UPLO = 'L') is
            exactly zero.
    < 0:    if the return value is -i, the i-th argument had an illegal value
   -100:    failed to allocate enough internal memory


  Further Details
  ===============

  The normal use for SKTRF is to compute the U T U^T or L T L^T
  decomposition of a skew-symmetric matrix with pivoting. This mode
  is chosen by setting MODE = 'N' ("normal" mode). The other
  use of SKTRF is the computation the Pfaffian of a skew-symmetric matrix,
  which only requires a partial computation of T, this mode is chosen
  by setting MODE = 'P' ("Pfaffian" mode).

  Normal mode (MODE = 'N'):
  ========================

  If UPLO = 'U', the U*T*U^T decomposition of A is computed. U is a
  upper triangular unit matrix with the additional constraint
  U(1:n-1,n) = 0, and T a tridiagonal matrix. The upper diagonal
  of T is stored on exit in A(i,i+1) for i = 1 .. n-1. The column
  U(1:i-1, i) is stored in A(1:i-1,i+1).

  If UPLO = 'L', the L*T*L^T decomposition of A is computed. L is a
  lower triangular unit matrix with the additional constraint
  L(2:n,1) = 0, and T a tridiagonal matrix. The lower diagonal
  of T is stored on exit in A(i+1,i) for i = 1 .. n-1. The column
  L(i+1:n, i) is stored in A(i+1:n,i-1).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  0   e   u2  u3  u4 )              (  0                  )
    (      0   e   u3  u4 )              (  e   0              )
    (          0   e   u4 )              (  l2  e   0          )
    (              0   e  )              (  l2  l3  e   0      )
    (                  0  )              (  l2  l3  l4  e   0  )

  where e denotes the off-diagonal elements of T, and ui (li)
  denotes an element of the i-th column of U (L).

  Pfaffian mode (MODE = 'P'):
  ==========================

  For computing the Pfaffian, it is enough to bring A into a partial
  tridiagonal form. In particular, assuming n even, it is enough to
  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
  only the off-diagonal entries in the odd columns (if UPLO = 'L')
  or in the even columns (if UPLU = 'U') are properly computed by SKTRF.

  If UPLO = 'U', the U*pT*U^T decomposition of A is computed. U is a
  upper triangular unit matrix with the additional constraint
  U(1:i-1,i) = 0 for even i, and pT a partially tridiagonal matrix.
  The entries in the odd rows of the upper diagonal of pT are stored
  on exit in A(i,i+1) for i odd. The column U(1:i-1, i) for odd i
  is stored in A(1:i-1,i+1).

  If UPLO = 'L', the L*pT*L^T decomposition of A is computed. L is a
  lower triangular unit matrix with the additional constraint
  L(i+1:n,i) = 0 for odd i, and pT a partially tridiagonal matrix.
  The entries in odd columns in the lower diagonal of pT are stored
  on exit in A(i+1,i) for i odd. The column L(i+1:n, i) for i odd
  is stored in A(i+1:n,i-1).

  The contents of A on exit are illustrated by the following examples
  with n = 6:

  if UPLO = 'U':                       if UPLO = 'L':

    (  0   e   x   u3  x   u5 )              (  0                    )
    (      0   x   u3  x   u5 )              (  e   0                )
    (          0   e   x   u5 )              (  l2  x   0            )
    (              0   x   u5 )              (  l2  x   e   0        )
    (                  0   e  )              (  l2  x   l4  x   0    )
    (                      0  )              (  l2  x   l4  x   e  0 )

  where e denotes the off-diagonal elements of T, ui (li)
  denotes an element of the i-th column of U (L), and x denotes an
  element not computed by SKTRF.

****************************************************************************/

#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

int sktrf_s(int N, float *A, int *IPIV,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;
  int info = 0;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(IPIV == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mode != 'N' && mode != 'P') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int lwork;
    float *work, qwork;

    /*workspace query*/
    lwork = -1;
    PFAPACK_ssktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   &qwork, &lwork, &info);

    if(info<0) printf("Haeh1");

    lwork = (int)qwork;

    work = (float *)malloc(sizeof(float) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (float *)malloc(sizeof(float) * lwork);
      if(!work) {
	return -100;
      }
    }

    PFAPACK_ssktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   work, &lwork, &info);

    if(info<0) printf("Haeh2 %d\n", info);

    free(work);
  }

  return info;
}

int sktrf_d(int N, double *A, int *IPIV,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;
  int info = 0;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(IPIV == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mode != 'N' && mode != 'P') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int lwork;
    double *work, qwork;

    /*workspace query*/
    lwork = -1;
    PFAPACK_dsktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   &qwork, &lwork, &info);

    if(info<0) printf("Haeh1");

    lwork = (int)qwork;

    work = (double *)malloc(sizeof(double) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (double *)malloc(sizeof(double) * lwork);
      if(!work) {
	return -100;
      }
    }

    PFAPACK_dsktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   work, &lwork, &info);

    if(info<0) printf("Haeh2 %d\n", info);

    free(work);
  }

  return info;
}

int sktrf_c(int N, floatcmplx *A, int *IPIV,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;
  int info = 0;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(IPIV == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mode != 'N' && mode != 'P') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int lwork;
    floatcmplx *work, qwork;

    /*workspace query*/
    lwork = -1;
    PFAPACK_csktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   &qwork, &lwork, &info);

    if(info<0) printf("Haeh1");

    lwork = (int)realf(qwork);

    work = (floatcmplx *)malloc(sizeof(floatcmplx) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (floatcmplx *)malloc(sizeof(floatcmplx) * lwork);
      if(!work) {
	return -100;
      }
    }

    PFAPACK_csktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   work, &lwork, &info);

    if(info<0) printf("Haeh2 %d\n", info);

    free(work);
  }

  return info;
}

int sktrf_z(int N, doublecmplx *A, int *IPIV,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;
  int info = 0;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(IPIV == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mode != 'N' && mode != 'P') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int lwork;
    doublecmplx *work, qwork;

    /*workspace query*/
    lwork = -1;
    PFAPACK_zsktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   &qwork, &lwork, &info);

    if(info<0) printf("Haeh1");

    lwork = (int)real(qwork);

    work = (doublecmplx *)malloc(sizeof(doublecmplx) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (doublecmplx *)malloc(sizeof(doublecmplx) * lwork);
      if(!work) {
	return -100;
      }
    }

    PFAPACK_zsktrf(UPLO, MODE, &N, A, &ldim, IPIV,
		   work, &lwork, &info);

    if(info<0) printf("Haeh2 %d\n", info);

    free(work);
  }

  return info;
}

#ifdef __cplusplus
}
#endif
