/**********************************************************************
  Purpose
  =======

  skbpf10_x computes the Pfaffian of a banded skew-symmetric matrix, taking
  special care to avoid numerical under- or overflow.
  (at the cost of possible additional round-off errors).

  Usage
  =====

  int skbpf10_x(int N, int KD, datatype *A, type *PFAFF,
                const char *UPLO)

  where
     datatype = float, double, floatcmplx, or doublecmplx
  for x = s, d, c, or z

  When using a C++ compiler, the trailing _x may be omitted.

  Arguments
  =========

  N       (input) int
          Size of the matrix A. N >= 0.

  KD      (input) int
          Number of super- (if UPLO = 'U') or subdiagonals (if UPLO =
          'L'). KD >= 0.

  A       (input/output) datatype *
          pointer to a memory block of size (KD+1)*N*sizeof(datatype)
          If N is odd, AB is not referenced.
          If N is even:
          On entry, A must contain a Fortran matrtix. The upper or
          lower triangle of the skew-symmetric band matrix A, stored
          in the first KD+1 rows of the array.  The j-th column of A
          is stored in the j-th column of the array AB as follows:
           if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
           if UPLO = 'L', AB(1+i-j,j) = A(i,j) for j<=i<=min(n,j+kd).
          (For a more convenient access, the functions/macros
           bandupper_fortran() and bandlower_fortran() are provided)
          On exit, AB is overwritten with values generated during the
          computation.

  PFAFF   (output) datatype *
          pointer to a an array datatype[2]
          The value of the Pfaffian in the form
          PFAFF[0] * (10^PFAFF[1]).

  UPLO    (input) char *
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

  Return value
  ============

  The return value of skbpfa_x indicates whether an error occured:
      0:    successful exit
    < 0:    if the return value is -i, the i-th argument had an illegal value
   -100:    failed to allocate enough internal memory


  Further Details
  ===============

  For odd-sized matrices, the Pfaffian is 0 by default, hence no
  computation is needed in this case. For even-sized matrices,
  the Pfaffian is computed by bringing the skew-symmetric matrix A into
  tridiagonal form T by a unitary congruence transformation:
  Q^H * A * Q^* = T.
  This transformation is computed by the routine SKBTRD (for further
  details see there)

****************************************************************************/

#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

int skbpf10_s(int N, int KD, float *A, float *PFAFF,
	    const char *UPLO)
{
  char uplo;

  uplo = toupper(UPLO[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(PFAFF == NULL) {
    return -4;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -5;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    float *work;

    work = (float *)malloc(sizeof(float) * (3*N-1));
    if(!work) {
      return -100;
    }

    PFAPACK_sskbpf10(UPLO, &N, &KD, A, &ldim, PFAFF,
		     work, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
  }
  else {
    PFAFF[0] = 1.0;
    PFAFF[1] = 0.0;
  }

  return 0;
}

int skbpf10_d(int N, int KD, double *A, double *PFAFF,
	    const char *UPLO)
{
  char uplo;

  uplo = toupper(UPLO[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(PFAFF == NULL) {
    return -4;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -5;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    double *work;

    work = (double *)malloc(sizeof(double) * (3*N-1));
    if(!work) {
      return -100;
    }

    PFAPACK_dskbpf10(UPLO, &N, &KD, A, &ldim, PFAFF,
		     work, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
  }
  else {
    PFAFF[0] = 1.0;
    PFAFF[1] = 0.0;
  }

  return 0;
}

int skbpf10_c(int N, int KD, floatcmplx *A, floatcmplx *PFAFF,
	      const char *UPLO)
{
  char uplo;

  uplo = toupper(UPLO[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(PFAFF == NULL) {
    return -4;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -5;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    floatcmplx *work;
    float *rwork;

    work = (floatcmplx *)malloc(sizeof(floatcmplx) * N);
    if(!work) {
      return -100;
    }

    rwork = (float *)malloc(sizeof(float) * (2*N-1));
    if(!rwork) {
      free(work);
      return -100;
    }

    PFAPACK_cskbpf10(UPLO, &N, &KD, A, &ldim, PFAFF,
		     work, rwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(rwork);
    free(work);
  }
  else {
    PFAFF[0] = floatcmplx_one;
    PFAFF[1] = floatcmplx_zero;
  }

  return 0;
}

int skbpf10_z(int N, int KD, doublecmplx *A, doublecmplx *PFAFF,
	      const char *UPLO)
{
  char uplo;

  uplo = toupper(UPLO[0]);

  if(N < 0) {
    return -1;
  }
  else if(KD < 0) {
    return -2;
  }
  else if(A == NULL) {
    return -3;
  }
  else if(PFAFF == NULL) {
    return -4;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -5;
  }

  if(N > 0) {
    int ldim = KD+1;
    int info = 0;
    doublecmplx *work;
    double *rwork;

    work = (doublecmplx *)malloc(sizeof(doublecmplx) * N);
    if(!work) {
      return -100;
    }

    rwork = (double *)malloc(sizeof(double) * (2*N-1));
    if(!rwork) {
      free(work);
      return -100;
    }

    PFAPACK_zskbpf10(UPLO, &N, &KD, A, &ldim, PFAFF,
		     work, rwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(rwork);
    free(work);
  }
  else {
    PFAFF[0] = doublecmplx_one;
    PFAFF[1] = doublecmplx_zero;
  }

  return 0;
}

#ifdef __cplusplus
}
#endif
