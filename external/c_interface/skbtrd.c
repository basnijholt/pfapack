/**********************************************************************
  Purpose
  =======

  skbtrd_x reduces a real or complex skew-symmetric band matrix A to real
  skew-symmetric tridiagonal form T by an unitary congruence
  transformation (which reduces to an orthogonal similarity
  transfromation in the real case):
     Q^dagger * A * Q^* = T.

  Usage
  =====

  int skbtrd_x(int N, int KD, float *A, float *DETQ, float *Q,
	     const char *UPLO, const char *MODE, const char *VECT)

  where
     datatype = float, double, floatcmplx, or doublecmplx
  for x = s, d, c, or z

  When using a C++ compiler, the trailing _x may be omitted.

  Arguments
  =========

  N       (input) int
          Size of the matrix A. N must be even, if MODE = 'P'. N >= 0.

  KD      (input) int
          Number of super- (if UPLO = 'U') or subdiagonals (if UPLO =
          'L'). KD >= 0.

  A       (input/output) datatype *,
          pointer to a memory block of size (KD+1)*N*sizeof(datatype)
          On entry, A must contain a Fortran array. The upper or
          lower triangle of the skew-symmetric band matrix A, stored
          in the first KD+1 rows of the array.  The j-th column of A
          is stored in the j-th column of the array AB as follows:
           if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
           if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd).
          (For a more convenient access, the functions/macros
           bandupper_fortran() and bandlower_fortran() are provided)
          On exit, the zero diagonal elements of AB are left unchanged,
          if KD > 0, the elements on the first superdiagonal (if UPLO =
          'U') or the first subdiagonal (if UPLO = 'L') are overwritten
          by the off-diagonal elements of T; the rest of AB is
          overwritten by values generated during the reduction. If
          MODE = 'P', only the off-diagonal entries in the odd rows
          (columns) are computed for UPLO = 'U' (UPLO = 'L').

  DETQ    (output) datatype *,
          pointer to variable of type datatype
          May be NULL, if datatype = float or double.
          The value of the determinant of Q, which is a
          pure phase factor (in the real case, DETQ=1 always).
          Computing DETQ does not require to form Q explicitly.

  Q       (input/output) datatype *,
          pointer to a memory block of size N*N*sizeof(datatype)
          Q must be a matrix in Fortran format. If VECT = 'N',
          it is not referenced and may be NULL.
          On entry:
            if VECT = 'U', then Q must contain an N-by-N matrix X;
            if VECT = 'V'  then Q need not be set.
          On exit:
             if VECT = 'V', Q contains the N-by-N unitary matrix Q;
             if VECT = 'U', Q contains the product X*Q;

  UPLO    (input) char *
             = 'U':  Upper triangle of A is stored;
             = 'L':  Lower triangle of A is stored.

  MODE    (input) char *
             = 'N':  A is fully tridiagonalized
             = 'P':  A is partially tridiagonalized for Pfaffian computation

  VECT    (input) char *
            = 'N':  do not form Q;
            = 'V':  form Q;
            = 'U':  update a matrix X, by forming X*Q.

  Return value
  ============

  The return value of sktrd_x indicates whether an error occured:
      0:    successful exit
    < 0:    if the return value is -i, the i-th argument had an illegal value
   -100:    failed to allocate enough internal memory


  Further Details
  ===============

  The storage scheme for the skew-symmetric matrix is identical to the
  storage scheme for symmetric or Hermitian band matrices in LAPACK,
  i.e. the diagonal and the KD super- or subdiagonals are stored in a
  Fortran array with KD+1 rows and N columns. Note that the zero
  diagonal must also be explicitely stored (this was done to keep the
  structure of the program identical to the symmetric case)

  In particular this means that if
  - UPLO = 'U', then AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j

    Example: N=5, KD=2

  (  0     a12   a13              )
  (  -a12  0     a23   a24        )
  (  -a13  -a23  0     a34   a35  )
  (        -a24  -a34  0     a45  )
  (              -a35  -a45  0    )

  is stored as

   x   x   a13 a24 a35
   x   a12 a23 a34 a45
   0   0   0   0   0

  where x denotes an unused entry

  - UPLO = 'L', then AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd)

    Example: N=5, KD=2

  (  0     -a21  -a31              )
  (  a21   0     -a32  -a42        )
  (  a31   a32   0     -a43  -a53  )
  (        a42   a43   0     -a54  )
  (              a53   a54   0     )

  is stored as

   0   0   0   0   0
   a21 a32 a43 a54 x
   a31 a42 a53 x   x

  where x denotes an unused entry

****************************************************************************/

#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

int skbtrd_s(int N, int KD, float *A, float *DETQ, float *Q,
	     const char *UPLO, const char *MODE, const char *VECT)
{
  char uplo, mode, vect;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);
  vect = toupper(VECT[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(Q == NULL && (vect == 'V' || vect == 'U')) {
    return -5;
  }
  else if(uplo != 'U' && uplo != 'L') {
    return -6;
  }
  else if(mode != 'N' && mode != 'P') {
    return -7;
  }
  else if(vect != 'N' && vect !='U' && vect != 'V') {
    return -8;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    float *work, *e;

    e = (float*)malloc(sizeof(float) * N);
    if(!e) {
      return -100;
    }

    work = (float *)malloc(sizeof(float) * 2*N);
    if(!work) {
      free(e);
      return -100;
    }

    PFAPACK_sskbtrd(VECT, UPLO, MODE, &N, &KD, A, &ldim, e, Q, &N,
		    work, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(e);
    free(work);
  }

  if(DETQ) *DETQ = 1;

  return 0;
}

int skbtrd_d(int N, int KD, double *A, double *DETQ, double *Q,
	     const char *UPLO, const char *MODE, const char *VECT)
{
  char uplo, mode, vect;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);
  vect = toupper(VECT[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(Q == NULL && (vect == 'V' || vect == 'U')) {
    return -5;
  }
  else if(uplo != 'U' && uplo != 'L') {
    return -6;
  }
  else if(mode != 'N' && mode != 'P') {
    return -7;
  }
  else if(vect != 'N' && vect !='U' && vect != 'V') {
    return -8;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    double *work, *e;

    e = (double*)malloc(sizeof(double) * N);
    if(!e) {
      return -100;
    }

    work = (double *)malloc(sizeof(double) * 2*N);
    if(!work) {
      free(e);
      return -100;
    }

    PFAPACK_dskbtrd(VECT, UPLO, MODE, &N, &KD, A, &ldim, e, Q, &N,
		    work, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(e);
    free(work);
  }

  if(DETQ) *DETQ = 1;

  return 0;
}

int skbtrd_c(int N, int KD, floatcmplx *A, floatcmplx *DETQ, floatcmplx *Q,
	     const char *UPLO, const char *MODE, const char *VECT)
{
  char uplo, mode, vect;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);
  vect = toupper(VECT[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(DETQ == NULL) {
    return -4;
  }
  else if(Q == NULL && (vect == 'V' || vect == 'U')) {
    return -5;
  }
  else if(uplo != 'U' && uplo != 'L') {
    return -6;
  }
  else if(mode != 'N' && mode != 'P') {
    return -7;
  }
  else if(vect != 'N' && vect !='U' && vect != 'V') {
    return -8;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    float *rwork, *e;
    floatcmplx *work;

    e = (float*)malloc(sizeof(float) * N);
    if(!e) {
      return -100;
    }

    rwork = (float *)malloc(sizeof(float) * N);
    if(!rwork) {
      free(e);
      return -100;
    }

    work = (floatcmplx *)malloc(sizeof(floatcmplx) * N);
    if(!work) {
      free(rwork);
      free(e);
      return -100;
    }

    PFAPACK_cskbtrd(VECT, UPLO, MODE, &N, &KD, A, &ldim, e, DETQ, Q, &N,
		    work, rwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
    free(rwork);
    free(e);
  }

  return 0;
}

int skbtrd_z(int N, int KD, doublecmplx *A, doublecmplx *DETQ, doublecmplx *Q,
	     const char *UPLO, const char *MODE, const char *VECT)
{
  char uplo, mode, vect;

  uplo = toupper(UPLO[0]);
  mode = toupper(MODE[0]);
  vect = toupper(VECT[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(DETQ == NULL) {
    return -4;
  }
  else if(Q == NULL && (vect == 'V' || vect == 'U')) {
    return -5;
  }
  else if(uplo != 'U' && uplo != 'L') {
    return -6;
  }
  else if(mode != 'N' && mode != 'P') {
    return -7;
  }
  else if(vect != 'N' && vect !='U' && vect != 'V') {
    return -8;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    double *rwork, *e;
    doublecmplx *work;

    e = (double*)malloc(sizeof(double) * N);
    if(!e) {
      return -100;
    }

    rwork = (double *)malloc(sizeof(double) * N);
    if(!rwork) {
      free(e);
      return -100;
    }

    work = (doublecmplx *)malloc(sizeof(doublecmplx) * N);
    if(!work) {
      free(rwork);
      free(e);
      return -100;
    }

    PFAPACK_zskbtrd(VECT, UPLO, MODE, &N, &KD, A, &ldim, e, DETQ, Q, &N,
		    work, rwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
    free(rwork);
    free(e);
  }

  return 0;
}

#ifdef __cplusplus
}
#endif
