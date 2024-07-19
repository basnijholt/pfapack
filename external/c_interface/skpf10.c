/**********************************************************************
  Purpose
  =======

  skpf10_x computes the Pfaffian of a dense skew-symmetric matrix, taking
  special care to avoid numerical under- or overflow.
  (at the cost of possible additional round-off errors).

  Usage
  =====

  int skpf10_x(int N, datatype *A, type *PFAFF, const char *UPLO,
               const char *MTHD)

  where
     datatype = float, double, floatcmplx, or doublecmplx
  for x = s, d, c, or z

  When using a C++ compiler, the trailing _x may be omitted.

  Arguments
  =========

  N       (input) int
          Size of the matrix A. N >= 0.

  A       (input/output) datatype *
          pointer to a memory block of size N*N*sizeof(datatype)
          On entry, the skew-symmetric matrix A.
             If UPLO = 'U', the upper triangular part of A contains
                the upper triangular part of the matrix A, and the
                strictly lower triangular part of A is not referenced.
             If UPLO = 'L', the lower triangular part of A contains
                the lower triangular part of the matrix A, and the
                strictly upper triangular part of A is not referenced.
          If the matrix size is odd, A is not referenced. If the matrix
          size is even, A is overwritten by values generated during
          the computation.

  PFAFF   (output) datatype *
          pointer to a an array datatype[2]
          The value of the Pfaffian in the form
          PFAFF[0] * (10^PFAFF[1]).

  UPLO    (input) char *
            = 'U':  Upper triangle of A is stored;
            = 'L':  Lower triangle of A is stored.

  MTHD    (input) char *
            = 'P': Compute Pfaffian using Parlett-Reid algorithm (recommended)
            = 'H': Compute Pfaffian using Householder reflections

  Return value
  ============

  The return value of skpfa_x indicates whether an error occured:
      0:    successful exit
    < 0:    if the return value is -i, the i-th argument had an illegal value
   -100:    failed to allocate enough internal memory


  Further Details
  ===============

  For odd-sized matrices, the Pfaffian is 0 by default, hence no
  computation is needed in this case. For even-sized matrices,
  the Pfaffian is computed by bringing the skew-symmetric matrix A into
  a partial tridiagonal form pT, either by computing a partial L pT L^T
  decomposition (MTHD = 'P'), or by a by a unitary congruence transformation
  Q^H * A * Q^* = pT (MTHD = 'H').
  These transformations are computed by the routines DSKTRF or DSKTRD,
  respectively (for further details see there).

****************************************************************************/


#include "commondefs.h"

#ifdef __cplusplus
extern "C" {
#endif

int skpf10_s(int N, float *A, float *PFAFF,
	    const char *UPLO, const char *MTHD)
{
  char uplo, mthd;

  uplo = toupper(UPLO[0]);
  mthd = toupper(MTHD[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(PFAFF == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mthd != 'P' && mthd != 'H') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int info = 0;
    int *iwork, lwork;
    float *work, qwork;

    iwork = (int *)malloc(sizeof(int) * N);
    if(!iwork) {
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_sskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, &qwork, &lwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)qwork;

    work = (float *)malloc(sizeof(float) * lwork);
    if(!work) {
      /*try minimal workspace*/
      if(mthd == 'P') lwork = 1;
      else lwork = 2*N-1;

      work = (float *)malloc(sizeof(float) * lwork);
      if(!work) {
	free(iwork);
	return -100;
      }
    }

    PFAPACK_sskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, work, &lwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
    free(iwork);
  }
  else {
    PFAFF[0] = 1.0;
    PFAFF[1] = 0.0;
  }

  return 0;
}

int skpf10_d(int N, double *A, double *PFAFF,
	    const char *UPLO, const char *MTHD)
{
  char uplo, mthd;

  uplo = toupper(UPLO[0]);
  mthd = toupper(MTHD[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(PFAFF == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mthd != 'P' && mthd != 'H') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int info = 0;
    int *iwork, lwork;
    double *work, qwork;

    iwork = (int *)malloc(sizeof(int) * N);
    if(!iwork) {
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_dskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, &qwork, &lwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)qwork;

    work = (double *)malloc(sizeof(double) * lwork);
    if(!work) {
      /*try minimal workspace*/
      if(mthd == 'P') lwork = 1;
      else lwork = 2*N-1;

      work = (double *)malloc(sizeof(double) * lwork);
      if(!work) {
	free(iwork);
	return -100;
      }
    }

    PFAPACK_dskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, work, &lwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
    free(iwork);
  }
  else {
    PFAFF[0] = 1.0;
    PFAFF[1] = 0.0;
  }

  return 0;
}

int skpf10_c(int N, floatcmplx *A, floatcmplx *PFAFF,
	     const char *UPLO, const char *MTHD)
{
  char uplo, mthd;

  uplo = toupper(UPLO[0]);
  mthd = toupper(MTHD[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(PFAFF == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mthd != 'P' && mthd != 'H') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int info = 0;
    int *iwork, lwork;
    float *rwork;
    floatcmplx *work, qwork;

    iwork = (int *)malloc(sizeof(int) * N);
    if(!iwork) {
      return -100;
    }
    rwork = (float *)malloc(sizeof(float) * (N-1));
    if(!rwork) {
      free(iwork);
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_cskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, &qwork, &lwork, rwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)realf(qwork);

    work = (floatcmplx *)malloc(sizeof(floatcmplx) * lwork);
    if(!work) {
      /*try minimal workspace*/
      if(mthd == 'P') lwork = 1;
      else lwork = 2*N-1;

      work = (floatcmplx *)malloc(sizeof(floatcmplx) * lwork);
      if(!work) {
	free(rwork);
	free(iwork);
	return -100;
      }
    }

    PFAPACK_cskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, work, &lwork, rwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
    free(rwork);
    free(iwork);
  }
  else {
    PFAFF[0] = floatcmplx_one;
    PFAFF[1] = floatcmplx_zero;
  }

  return 0;
}

int skpf10_z(int N, doublecmplx *A, doublecmplx *PFAFF,
	     const char *UPLO, const char *MTHD)
{
  char uplo, mthd;

  uplo = toupper(UPLO[0]);
  mthd = toupper(MTHD[0]);

  if(N<0) {
    return -1;
  }
  else if(A == NULL) {
    return -2;
  }
  else if(PFAFF == NULL) {
    return -3;
  }
  else if(uplo != 'U' && uplo !='L') {
    return -4;
  }
  else if(mthd != 'P' && mthd != 'H') {
    return -5;
  }

  if(N > 0) {
    int ldim = N;
    int info = 0;
    int *iwork, lwork;
    double *rwork;
    doublecmplx *work, qwork;

    iwork = (int *)malloc(sizeof(int) * N);
    if(!iwork) {
      return -100;
    }
    rwork = (double *)malloc(sizeof(double) * (N-1));
    if(!rwork) {
      free(iwork);
      return -100;
    }

    /*workspace query*/
    lwork = -1;
    PFAPACK_zskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, &qwork, &lwork, rwork, &info);

    if(info) printf("Haeh1");

    lwork = (int)real(qwork);

    work = (doublecmplx *)malloc(sizeof(doublecmplx) * lwork);
    if(!work) {
      /*try minimal workspace*/
      if(mthd == 'P') lwork = 1;
      else lwork = 2*N-1;

      work = (doublecmplx *)malloc(sizeof(doublecmplx) * lwork);
      if(!work) {
	free(rwork);
	free(iwork);
	return -100;
      }
    }

    PFAPACK_zskpf10(UPLO, MTHD, &N, A, &ldim, PFAFF,
		    iwork, work, &lwork, rwork, &info);

    if(info) printf("Haeh2 %d\n", info);

    free(work);
    free(rwork);
    free(iwork);
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
