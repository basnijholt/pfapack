/**********************************************************************
  Purpose
  =======

  sktrd_x reduces a real or complex skew-symmetric matrix A to real
  skew-symmetric tridiagonal form T by a unitary congruence
  transformation (which reduces to an orthogonal similarity
  transformation in the real case):
  Q^T * A * Q = T.

  This subroutine uses, if enough memory is available, the
  blocked version of the Householder algorithm.

  Usage
  =====

  int sktrd_x(int N, datatype *A, datatype *TAU,
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
          On exit, if MODE = 'N':
            If UPLO = 'U', the diagonal and first superdiagonal
              of A are overwritten by the corresponding elements of the
              tridiagonal matrix T, and the elements above the first
              superdiagonal, with the array TAU, represent the unitary
              matrix Q as a product of elementary reflectors;
            If UPLO = 'L', the diagonal and first subdiagonal of A are over-
              written by the corresponding elements of the tridiagonal
              matrix T, and the elements below the first subdiagonal, with
              the array TAU, represent the unitary matrix Q as a product
              of elementary reflectors.
            See Further Details, also for information about MODE = 'P'.

  TAU     (output) datatype *
          pointer to array of dimesniuon N-1
          The scalar factors of the elementary reflectors (see Further
          Details).

  UPLO    (input) char *
             = 'U':  Upper triangle of A is stored;
             = 'L':  Lower triangle of A is stored.

  MODE    (input) char *
             = 'N':  A is fully tridiagonalized
             = 'P':  A is partially tridiagonalized for Pfaffian computation

  Return value
  ============

  The return value of sktrd_x indicates whether an error occured:
      0:    successful exit
    < 0:    if the return value is -i, the i-th argument had an illegal value
   -100:    failed to allocate enough internal memory

  Further Details
  ===============

  The normal use for SKTRD is to compute the tridiagonal form of
  a skew-symmetric matrix under an orthogonal similarity transformation,
  and chosen by setting MODE = 'N' ("normal" mode). The other
  use of SKTRD is the computation the Pfaffian of a skew-symmetric matrix,
  which only requires a partial tridiagonalization, this mode is chosen
  by setting MODE = 'P' ("Pfaffian" mode).

  Normal mode (MODE = 'N'):
  ========================

  The routine computes a tridiagonal matrix T and an orthogonal Q such
  that A = Q * T * Q^T .

  If UPLO = 'U', the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-1) . . . H(2) H(1).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  A(1:i-1,i+1), and tau in TAU(i).

  If UPLO = 'L', the matrix Q is represented as a product of elementary
  reflectors

     Q = H(1) H(2) . . . H(n-1).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  0   e   v2  v3  v4 )              (  0                  )
    (      0   e   v3  v4 )              (  e   0              )
    (          0   e   v4 )              (  v1  e   0          )
    (              0   e  )              (  v1  v2  e   0      )
    (                  0  )              (  v1  v2  v3  e   0  )

  where d and e denote diagonal and off-diagonal elements of T, and vi
  denotes an element of the vector defining H(i).

  The LAPACK routine DORGTR can be used to form the transformation
  matrix explicitely, and DORMTR can be used to multiply another
  matrix without forming the transformation.

  Pfaffian mode (MODE = 'P'):
  ==========================

  For computing the Pfaffian, it is enough to bring A into a partial
  tridiagonal form. In particular, assuming n even, it is enough to
  bring A into a form with A(i,j) = A(j,i) = 0 for i > j+1 with j odd
  (this is computed if UPLO = 'L'), or A(i,j) = A(j,i) = 0 for
  i > j-1 with j even (this is computed if UPLO = 'U'). Note that
  only the off-diagonal entries in the odd columns (if UPLO = 'L')
  or in the even columns (if UPLU = 'U') are properly computed by SKTRD.

  A is brought into this special form pT using an orthogonal matrix Q:
  A = Q * pT * Q^T

  If UPLO = 'U', the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-1) H(n-3) . . . H(3) H(1).

  Each H(i) has the form

     H(i) = I - tau * v * v^T

  where tau is a real scalar, and v is a real vector with
  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  A(1:i-1,i+1), and tau in TAU(i).

  If UPLO = 'L', the matrix Q is represented as a product of elementary
  reflectors

     Q = H(1) H(3) . . . H(n-3) H(n-1).

  Each H(i) has the form

     H(i) = I - tau * v * v^T

  where tau is a real scalar, and v is a real vector with
  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 6:

  if UPLO = 'U':                       if UPLO = 'L':

    (  0   e   x   v3  x   v5 )        (  0                      )
    (      0   x   v3  x   v5 )        (  e   0                  )
    (          0   e   x   v5 )        (  v1  x   0              )
    (              0   x   v5 )        (  v1  x   e   0          )
    (                  0   e  )        (  v1  x   v3  x   0      )
    (                      0  )        (  v1  x   v3  x   e   0  )

  where d and e denote diagonal and off-diagonal elements of T, vi
  denotes an element of the vector defining H(i), and x denotes an
  element not computed by SKTRD.

****************************************************************************/

#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

int sktrd_s(int N, float *A, float *TAU,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(TAU == NULL) {
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
    int info = 0;
    int lwork;
    float *work, *e, qwork;

    e = (float *)malloc(sizeof(float) * (N-1));
    if(!e) {
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_ssktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   &qwork, &lwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)qwork;

    work = (float *)malloc(sizeof(float) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (float *)malloc(sizeof(float) * lwork);
      if(!work) {
	free(e);
	return -100;
      }
    }

    PFAPACK_ssktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   work, &lwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(e);
    free(work);
  }

  return 0;
}

int sktrd_d(int N, double *A, double *TAU,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(TAU == NULL) {
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
    int info = 0;
    int lwork;
    double *work, *e, qwork;

    e = (double *)malloc(sizeof(double) * (N-1));
    if(!e) {
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_dsktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   &qwork, &lwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)qwork;

    work = (double *)malloc(sizeof(double) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (double *)malloc(sizeof(double) * lwork);
      if(!work) {
	free(e);
	return -100;
      }
    }

    PFAPACK_dsktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   work, &lwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(e);
    free(work);
  }

  return 0;
}

int sktrd_c(int N, floatcmplx *A, floatcmplx *TAU,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(TAU == NULL) {
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
    int info = 0;
    int lwork;
    float *e;
    floatcmplx *work, qwork;

    e = (float *)malloc(sizeof(float) * (N-1));
    if(!e) {
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_csktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   &qwork, &lwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)realf(qwork);

    work = (floatcmplx *)malloc(sizeof(floatcmplx) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (floatcmplx *)malloc(sizeof(floatcmplx) * lwork);
      if(!work) {
	free(e);
	return -100;
      }
    }

    PFAPACK_csktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   work, &lwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(e);
    free(work);
  }

  return 0;
}

int sktrd_z(int N, doublecmplx *A, doublecmplx *TAU,
	    const char *UPLO, const char *MODE)
{
  char uplo, mode;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(TAU == NULL) {
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
    int info = 0;
    int lwork;
    double *e;
    doublecmplx *work, qwork;

    e = (double *)malloc(sizeof(double) * (N-1));
    if(!e) {
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_zsktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   &qwork, &lwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)real(qwork);

    work = (doublecmplx *)malloc(sizeof(doublecmplx) * lwork);
    if(!work) {
      /*try minimal workspace*/
      lwork = 1;

      work = (doublecmplx *)malloc(sizeof(doublecmplx) * lwork);
      if(!work) {
	free(e);
	return -100;
      }
    }

    PFAPACK_zsktrd(UPLO, MODE, &N, A, &ldim, e, TAU,
		   work, &lwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(e);
    free(work);
  }

  return 0;
}

#ifdef __cplusplus
}
#endif
