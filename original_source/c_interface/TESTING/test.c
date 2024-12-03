#include "pfapack.h"

#ifdef __cplusplus

#include <cstdio>
#include <cstdlib>
#include <cstring>

using std::printf;
using std::rand;
using std::malloc;
using std::free;
using std::memset;
using std::memcpy;

static const floatcmplx floatcmplx_one = std::complex<float>(1);
static const doublecmplx doublecmplx_one = std::complex<double>(1);

static const floatcmplx floatcmplx_zero = std::complex<float>(0);
static const doublecmplx doublecmplx_zero = std::complex<double>(0);

#else

#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#if __STDC_VERSION__ >= 199901L

static const floatcmplx floatcmplx_one = 1;
static const doublecmplx doublecmplx_one = 1;

static const floatcmplx floatcmplx_zero = 0;
static const doublecmplx doublecmplx_zero = 0;

#else

static const floatcmplx floatcmplx_one = { 1, 0 };
static const doublecmplx doublecmplx_one = { 1, 0 };

static const floatcmplx floatcmplx_zero = { 0, 0 };
static const doublecmplx doublecmplx_zero = { 0, 0 };

#endif

#endif

void make_skew_mat_s(int N, float *A, double prob)
{
  int i, j;

  memset(A, 0, sizeof(float)*N*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if( rand()/(double)RAND_MAX < prob) {
	A[dense_fortran(i, j, N)] = rand()*2.0/RAND_MAX-1.0;
	A[dense_fortran(j, i, N)] = -A[dense_fortran(i, j, N)];
      }
    }
  }
}

void make_skew_mat_d(int N, double *A, double prob)
{
  int i, j;

  memset(A, 0, sizeof(double)*N*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if( rand()/(double)RAND_MAX < prob) {
	A[dense_fortran(i, j, N)] = rand()*2.0/RAND_MAX-1.0;
	A[dense_fortran(j, i, N)] = -A[dense_fortran(i, j, N)];
      }
    }
  }
}

void make_skew_mat_c(int N, floatcmplx *A, double prob)
{
  int i, j;

  memset(A, 0, sizeof(floatcmplx)*N*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if( rand()/(double)RAND_MAX < prob) {
#ifdef __cplusplus
	A[dense_fortran(i, j, N)] = floatcmplx(rand()*2.0/RAND_MAX-1.0,
					       rand()*2.0/RAND_MAX-1.0);
	A[dense_fortran(j, i, N)] = -A[dense_fortran(i, j, N)];
#elif __STDC_VERSION__ >= 199901L
	A[dense_fortran(i, j, N)] = (rand()*2.0/RAND_MAX-1.0) +
	  I*(rand()*2.0/RAND_MAX-1.0);
	A[dense_fortran(j, i, N)] = -A[dense_fortran(i, j, N)];
#else
	A[dense_fortran(i, j, N)].re = rand()*2.0/RAND_MAX-1.0;
	A[dense_fortran(i, j, N)].im = rand()*2.0/RAND_MAX-1.0;
	A[dense_fortran(j, i, N)].re = -A[dense_fortran(i, j, N)].re;
	A[dense_fortran(j, i, N)].im = -A[dense_fortran(i, j, N)].im;
#endif
      }
    }
  }
}

void make_skew_mat_z(int N, doublecmplx *A, double prob)
{
  int i, j;

  memset(A, 0, sizeof(doublecmplx)*N*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if( rand()/(double)RAND_MAX < prob) {
#ifdef __cplusplus
	A[dense_fortran(i, j, N)] = doublecmplx(rand()*2.0/RAND_MAX-1.0,
						rand()*2.0/RAND_MAX-1.0);
	A[dense_fortran(j, i, N)] = -A[dense_fortran(i, j, N)];
#elif __STDC_VERSION__ >= 199901L
	A[dense_fortran(i, j, N)] = (rand()*2.0/RAND_MAX-1.0) +
	  I*(rand()*2.0/RAND_MAX-1.0);
	A[dense_fortran(j, i, N)] = -A[dense_fortran(i, j, N)];
#else
	A[dense_fortran(i, j, N)].re = rand()*2.0/RAND_MAX-1.0;
	A[dense_fortran(i, j, N)].im = rand()*2.0/RAND_MAX-1.0;
	A[dense_fortran(j, i, N)].re = -A[dense_fortran(i, j, N)].re;
	A[dense_fortran(j, i, N)].im = -A[dense_fortran(i, j, N)].im;
#endif
      }
    }
  }
}

/*make random baned matrices*/
void make_skew_mat_banded_s(int N, int KD, float *Au,
			    float *Ad, double prob)
{
  int i,j;

  if(Au) memset(Au, 0, sizeof(float)*(KD+1)*N);
  if(Ad) memset(Ad, 0, sizeof(float)*(KD+1)*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if(i-j < KD) {
	if( rand()/(double)RAND_MAX < prob) {
	  float tmp = rand()*2.0/RAND_MAX-1.0;

	  if(Au) Au[bandupper_fortran(j, i, N, KD)] = tmp;
	  if(Ad) Ad[bandlower_fortran(i, j, N, KD)] = -tmp;
	}
      }
    }
  }
}

void make_skew_mat_banded_d(int N, int KD, double *Au,
			    double *Ad, double prob)
{
  int i,j;

  if(Au) memset(Au, 0, sizeof(double)*(KD+1)*N);
  if(Ad) memset(Ad, 0, sizeof(double)*(KD+1)*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if(i-j < KD) {
	if( rand()/(double)RAND_MAX < prob) {
	  double tmp = rand()*2.0/RAND_MAX-1.0;

	  if(Au) Au[bandupper_fortran(j, i, N, KD)] = tmp;
	  if(Ad) Ad[bandlower_fortran(i, j, N, KD)] = -tmp;
	}
      }
    }
  }
}

void make_skew_mat_banded_c(int N, int KD, floatcmplx *Au,
			    floatcmplx *Ad, double prob)
{
  int i,j;

  if(Au) memset(Au, 0, sizeof(floatcmplx)*(KD+1)*N);
  if(Ad) memset(Ad, 0, sizeof(floatcmplx)*(KD+1)*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if(i-j < KD) {
	if( rand()/(double)RAND_MAX < prob) {
	  float tmpre = rand()*2.0/RAND_MAX-1.0;
	  float tmpim = rand()*2.0/RAND_MAX-1.0;

	  if(Au) {
            #ifdef __cplusplus
	    Au[bandupper_fortran(j, i, N, KD)] = floatcmplx(tmpre,
							    tmpim);
            #elif __STDC_VERSION__ >= 199901L
	    Au[bandupper_fortran(j, i, N, KD)] = tmpre + I*tmpim;
            #else
	    Au[bandupper_fortran(j, i, N, KD)].re = tmpre;
	    Au[bandupper_fortran(j, i, N, KD)].im = tmpim;
	    #endif
	  }
	  if(Ad) {
	    #ifdef __cplusplus
	    Ad[bandlower_fortran(i, j, N, KD)] = -floatcmplx(tmpre,
							     tmpim);
            #elif __STDC_VERSION__ >= 199901L
	    Ad[bandlower_fortran(i, j, N, KD)] = -tmpre - I*tmpim;
            #else
	    Ad[bandlower_fortran(i, j, N, KD)].re = -tmpre;
	    Ad[bandlower_fortran(i, j, N, KD)].im = -tmpim;
	    #endif
	  }
	}
      }
    }
  }
}

void make_skew_mat_banded_z(int N, int KD, doublecmplx *Au,
			    doublecmplx *Ad, double prob)
{
  int i,j;

  if(Au) memset(Au, 0, sizeof(doublecmplx)*(KD+1)*N);
  if(Ad) memset(Ad, 0, sizeof(doublecmplx)*(KD+1)*N);

  for(i=1; i<=N; i++) {
    for(j=1; j<=i-1; j++) {
      if(i-j < KD) {
	if( rand()/(double)RAND_MAX < prob) {
	  double tmpre = rand()*2.0/RAND_MAX-1.0;
	  double tmpim = rand()*2.0/RAND_MAX-1.0;

	  if(Au) {
            #ifdef __cplusplus
	    Au[bandupper_fortran(j, i, N, KD)] = doublecmplx(tmpre,
							    tmpim);
            #elif __STDC_VERSION__ >= 199901L
	    Au[bandupper_fortran(j, i, N, KD)] = tmpre + I*tmpim;
            #else
	    Au[bandupper_fortran(j, i, N, KD)].re = tmpre;
	    Au[bandupper_fortran(j, i, N, KD)].im = tmpim;
	    #endif
	  }
	  if(Ad) {
	    #ifdef __cplusplus
	    Ad[bandlower_fortran(i, j, N, KD)] = -doublecmplx(tmpre,
							     tmpim);
            #elif __STDC_VERSION__ >= 199901L
	    Ad[bandlower_fortran(i, j, N, KD)] = -tmpre - I*tmpim;
            #else
	    Ad[bandlower_fortran(i, j, N, KD)].re = -tmpre;
	    Ad[bandlower_fortran(i, j, N, KD)].im = -tmpim;
	    #endif
	  }
	}
      }
    }
  }
}

#include "check.h"

#define MAX_SIZE 20

void check_pfaff_s()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      float *A, *Acpy;
      float pfaff, pfaff10[2];

      A = (float*)malloc(sizeof(float)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (float*)malloc(sizeof(float)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }

      /*check skpfa_x for various parameters*/
      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa_s(N, A, &pfaff, "U", "P");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa_s(N, A, &pfaff, "L", "P");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa_s(N, A, &pfaff, "U", "H");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa_s(N, A, &pfaff, "L", "H");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa(N, A, &pfaff, "U", "P");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa(N, A, &pfaff, "L", "P");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa(N, A, &pfaff, "U", "H");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpfa(N, A, &pfaff, "L", "H");
      F95_check_skpfa_s(&N, Acpy, A, &pfaff, "L", "H");

      #endif

      /*check skpf10_x for various parameters*/
      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10_s(N, A, pfaff10, "U", "P");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10_s(N, A, pfaff10, "L", "P");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10_s(N, A, pfaff10, "U", "H");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10_s(N, A, pfaff10, "L", "H");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10(N, A, pfaff10, "U", "P");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10(N, A, pfaff10, "L", "P");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10(N, A, pfaff10, "U", "H");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      skpf10(N, A, pfaff10, "L", "H");
      F95_check_skpf10_s(&N, Acpy, A, pfaff10, "L", "H");

      #endif

      free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	float *A, *Acpy;
	float pfaff, pfaff10[2];

	A = (float*)malloc(sizeof(float)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (float*)malloc(sizeof(float)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpfa_s(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_s(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpfa_s(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_s(&N, &KD, Acpy, A, &pfaff, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_s(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_s(&N, &KD, Acpy, A, &pfaff, "L");

	#endif

	/*check skbpf10_x for various parameters*/
	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpf10_s(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_s(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpf10_s(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_s(&N, &KD, Acpy, A, pfaff10, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_s(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_s(&N, &KD, Acpy, A, pfaff10, "L");

	#endif

	free(Acpy); free(A);
      }
    }
  }
}

void check_pfaff_d()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      double *A, *Acpy;
      double pfaff, pfaff10[2];

      A = (double*)malloc(sizeof(double)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (double*)malloc(sizeof(double)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }

      /*check skpfa_x for various parameters*/
      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa_d(N, A, &pfaff, "U", "P");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa_d(N, A, &pfaff, "L", "P");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa_d(N, A, &pfaff, "U", "H");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa_d(N, A, &pfaff, "L", "H");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa(N, A, &pfaff, "U", "P");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa(N, A, &pfaff, "L", "P");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa(N, A, &pfaff, "U", "H");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpfa(N, A, &pfaff, "L", "H");
      F95_check_skpfa_d(&N, Acpy, A, &pfaff, "L", "H");

      #endif

      /*check skpf10_x for various parameters*/
      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10_d(N, A, pfaff10, "U", "P");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10_d(N, A, pfaff10, "L", "P");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10_d(N, A, pfaff10, "U", "H");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10_d(N, A, pfaff10, "L", "H");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10(N, A, pfaff10, "U", "P");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10(N, A, pfaff10, "L", "P");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10(N, A, pfaff10, "U", "H");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      skpf10(N, A, pfaff10, "L", "H");
      F95_check_skpf10_d(&N, Acpy, A, pfaff10, "L", "H");

      #endif

      free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	double *A, *Acpy;
	double pfaff, pfaff10[2];

	A = (double*)malloc(sizeof(double)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (double*)malloc(sizeof(double)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpfa_d(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_d(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpfa_d(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_d(&N, &KD, Acpy, A, &pfaff, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_d(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_d(&N, &KD, Acpy, A, &pfaff, "L");

	#endif

	/*check skbpf10_x for various parameters*/
	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpf10_d(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_d(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpf10_d(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_d(&N, &KD, Acpy, A, pfaff10, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_d(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_d(&N, &KD, Acpy, A, pfaff10, "L");

	#endif

	free(Acpy); free(A);
      }
    }
  }
}

void check_pfaff_c()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      floatcmplx *A, *Acpy;
      floatcmplx pfaff, pfaff10[2];

      A = (floatcmplx*)malloc(sizeof(floatcmplx)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (floatcmplx*)malloc(sizeof(floatcmplx)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }

      /*check skpfa_x for various parameters*/
      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa_c(N, A, &pfaff, "U", "P");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa_c(N, A, &pfaff, "L", "P");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa_c(N, A, &pfaff, "U", "H");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa_c(N, A, &pfaff, "L", "H");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa(N, A, &pfaff, "U", "P");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa(N, A, &pfaff, "L", "P");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa(N, A, &pfaff, "U", "H");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpfa(N, A, &pfaff, "L", "H");
      F95_check_skpfa_c(&N, Acpy, A, &pfaff, "L", "H");

      #endif

      /*check skpf10_x for various parameters*/
      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10_c(N, A, pfaff10, "U", "P");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10_c(N, A, pfaff10, "L", "P");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10_c(N, A, pfaff10, "U", "H");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10_c(N, A, pfaff10, "L", "H");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10(N, A, pfaff10, "U", "P");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10(N, A, pfaff10, "L", "P");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10(N, A, pfaff10, "U", "H");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      skpf10(N, A, pfaff10, "L", "H");
      F95_check_skpf10_c(&N, Acpy, A, pfaff10, "L", "H");

      #endif

      free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	floatcmplx *A, *Acpy;
	floatcmplx pfaff, pfaff10[2];

	A = (floatcmplx*)malloc(sizeof(floatcmplx)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (floatcmplx*)malloc(sizeof(floatcmplx)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpfa_c(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_c(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpfa_c(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_c(&N, &KD, Acpy, A, &pfaff, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_c(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_c(&N, &KD, Acpy, A, &pfaff, "L");

	#endif

	/*check skbpf10_x for various parameters*/
	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpf10_c(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_c(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpf10_c(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_c(&N, &KD, Acpy, A, pfaff10, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_c(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_c(&N, &KD, Acpy, A, pfaff10, "L");

	#endif

	free(Acpy); free(A);
      }
    }
  }
}

void check_pfaff_z()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      doublecmplx *A, *Acpy;
      doublecmplx pfaff, pfaff10[2];

      A = (doublecmplx*)malloc(sizeof(doublecmplx)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (doublecmplx*)malloc(sizeof(doublecmplx)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }

      /*check skpfa_x for various parameters*/
      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa_z(N, A, &pfaff, "U", "P");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa_z(N, A, &pfaff, "L", "P");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa_z(N, A, &pfaff, "U", "H");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa_z(N, A, &pfaff, "L", "H");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa(N, A, &pfaff, "U", "P");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "U", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa(N, A, &pfaff, "L", "P");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "L", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa(N, A, &pfaff, "U", "H");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "U", "H");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpfa(N, A, &pfaff, "L", "H");
      F95_check_skpfa_z(&N, Acpy, A, &pfaff, "L", "H");

      #endif

      /*check skpf10_x for various parameters*/
      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10_z(N, A, pfaff10, "U", "P");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10_z(N, A, pfaff10, "L", "P");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10_z(N, A, pfaff10, "U", "H");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10_z(N, A, pfaff10, "L", "H");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "L", "H");

      #ifdef __cplusplus

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10(N, A, pfaff10, "U", "P");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "U", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10(N, A, pfaff10, "L", "P");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "L", "P");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10(N, A, pfaff10, "U", "H");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "U", "H");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      skpf10(N, A, pfaff10, "L", "H");
      F95_check_skpf10_z(&N, Acpy, A, pfaff10, "L", "H");

      #endif

      free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	doublecmplx *A, *Acpy;
	doublecmplx pfaff, pfaff10[2];

	A = (doublecmplx*)malloc(sizeof(doublecmplx)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (doublecmplx*)malloc(sizeof(doublecmplx)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpfa_z(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_z(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpfa_z(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_z(&N, &KD, Acpy, A, &pfaff, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "U");
	F95_check_skbpfa_z(&N, &KD, Acpy, A, &pfaff, "U");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpfa(N, KD, A, &pfaff, "L");
	F95_check_skbpfa_z(&N, &KD, Acpy, A, &pfaff, "L");

	#endif

	/*check skbpf10_x for various parameters*/
	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpf10_z(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_z(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpf10_z(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_z(&N, &KD, Acpy, A, pfaff10, "L");

	#ifdef __cplusplus

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "U");
	F95_check_skbpf10_z(&N, &KD, Acpy, A, pfaff10, "U");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbpf10(N, KD, A, pfaff10, "L");
	F95_check_skbpf10_z(&N, &KD, Acpy, A, pfaff10, "L");

	#endif

	free(Acpy); free(A);
      }
    }
  }
}

void check_decomp_s()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      float *A, *Acpy;
      float *tau;
      int *ipiv;

      A = (float*)malloc(sizeof(float)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (float*)malloc(sizeof(float)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }
      tau = (float*)malloc(sizeof(float)*(N-1));
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(Acpy); free(A);
	exit(10);
      }
      ipiv = (int *)malloc(sizeof(int)*N);
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(tau), free(Acpy); free(A);
	exit(10);
      }

      /*check sktrf_x for various parameters*/
      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrf_s(N, A, ipiv, "U", "N");
      F95_check_sktrf_s(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrf_s(N, A, ipiv, "L", "N");
      F95_check_sktrf_s(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrf_s(N, A, ipiv, "U", "P");
	F95_check_sktrf_s(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrf_s(N, A, ipiv, "L", "P");
	F95_check_sktrf_s(&N, Acpy, A, ipiv, "L", "P");
      }

      #ifdef __cplusplus

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrf(N, A, ipiv, "U", "N");
      F95_check_sktrf_s(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrf(N, A, ipiv, "L", "N");
      F95_check_sktrf_s(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrf(N, A, ipiv, "U", "P");
	F95_check_sktrf_s(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrf(N, A, ipiv, "L", "P");
	F95_check_sktrf_s(&N, Acpy, A, ipiv, "L", "P");
      }

      #endif

      /*check sktrd_x for various parameters*/
      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrd_s(N, A, tau, "U", "N");
      F95_check_sktrd_s(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrd_s(N, A, tau, "L", "N");
      F95_check_sktrd_s(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrd_s(N, A, tau, "U", "P");
	F95_check_sktrd_s(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrd_s(N, A, tau, "L", "P");
	F95_check_sktrd_s(&N, Acpy, A, tau, "L", "P");
      }

      #ifdef __cplusplus
      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrd(N, A, tau, "U", "N");
      F95_check_sktrd_s(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_s(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(float)*N*N);

      sktrd(N, A, tau, "L", "N");
      F95_check_sktrd_s(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrd(N, A, tau, "U", "P");
	F95_check_sktrd_s(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_s(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*N);

	sktrd(N, A, tau, "L", "P");
	F95_check_sktrd_s(&N, Acpy, A, tau, "L", "P");
      }

      #endif

      free(ipiv); free(tau); free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	float *A, *Acpy;
	float detQ, *Q;
	int i;

	A = (float*)malloc(sizeof(float)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (float*)malloc(sizeof(float)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}
	Q = (float*)malloc(sizeof(float)*N*N);
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(Acpy); free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd_s(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd_s(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd_s(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd_s(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd_s(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd_s(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd_s(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd_s(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd_s(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd_s(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd_s(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd_s(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#ifdef __cplusplus
	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_s(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_s(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(float)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_s(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#endif

	free(Q); free(Acpy); free(A);
      }
    }
  }
}

void check_decomp_d()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      double *A, *Acpy;
      double *tau;
      int *ipiv;

      A = (double*)malloc(sizeof(double)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (double*)malloc(sizeof(double)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }
      tau = (double*)malloc(sizeof(double)*(N-1));
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(Acpy); free(A);
	exit(10);
      }
      ipiv = (int *)malloc(sizeof(int)*N);
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(tau), free(Acpy); free(A);
	exit(10);
      }

      /*check sktrf_x for various parameters*/
      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrf_d(N, A, ipiv, "U", "N");
      F95_check_sktrf_d(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrf_d(N, A, ipiv, "L", "N");
      F95_check_sktrf_d(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrf_d(N, A, ipiv, "U", "P");
	F95_check_sktrf_d(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrf_d(N, A, ipiv, "L", "P");
	F95_check_sktrf_d(&N, Acpy, A, ipiv, "L", "P");
      }

      #ifdef __cplusplus

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrf(N, A, ipiv, "U", "N");
      F95_check_sktrf_d(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrf(N, A, ipiv, "L", "N");
      F95_check_sktrf_d(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrf(N, A, ipiv, "U", "P");
	F95_check_sktrf_d(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrf(N, A, ipiv, "L", "P");
	F95_check_sktrf_d(&N, Acpy, A, ipiv, "L", "P");
      }

      #endif

      /*check sktrd_x for various parameters*/
      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrd_d(N, A, tau, "U", "N");
      F95_check_sktrd_d(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrd_d(N, A, tau, "L", "N");
      F95_check_sktrd_d(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrd_d(N, A, tau, "U", "P");
	F95_check_sktrd_d(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrd_d(N, A, tau, "L", "P");
	F95_check_sktrd_d(&N, Acpy, A, tau, "L", "P");
      }

      #ifdef __cplusplus
      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrd(N, A, tau, "U", "N");
      F95_check_sktrd_d(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_d(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(double)*N*N);

      sktrd(N, A, tau, "L", "N");
      F95_check_sktrd_d(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrd(N, A, tau, "U", "P");
	F95_check_sktrd_d(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_d(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*N);

	sktrd(N, A, tau, "L", "P");
	F95_check_sktrd_d(&N, Acpy, A, tau, "L", "P");
      }

      #endif

      free(ipiv); free(tau); free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	double *A, *Acpy;
	double detQ, *Q;
	int i;

	A = (double*)malloc(sizeof(double)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (double*)malloc(sizeof(double)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}
	Q = (double*)malloc(sizeof(double)*N*N);
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(Acpy); free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd_d(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd_d(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd_d(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd_d(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd_d(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd_d(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd_d(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd_d(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd_d(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd_d(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd_d(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd_d(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#ifdef __cplusplus
	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = 1.0;

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_d(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_d(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(double)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = 1.0;

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_d(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#endif

	free(Q); free(Acpy); free(A);
      }
    }
  }
}

void check_decomp_c()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      floatcmplx *A, *Acpy;
      floatcmplx *tau;
      int *ipiv;

      A = (floatcmplx*)malloc(sizeof(floatcmplx)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (floatcmplx*)malloc(sizeof(floatcmplx)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }
      tau = (floatcmplx*)malloc(sizeof(floatcmplx)*(N-1));
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(Acpy); free(A);
	exit(10);
      }
      ipiv = (int *)malloc(sizeof(int)*N);
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(tau), free(Acpy); free(A);
	exit(10);
      }

      /*check sktrf_x for various parameters*/
      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrf_c(N, A, ipiv, "U", "N");
      F95_check_sktrf_c(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrf_c(N, A, ipiv, "L", "N");
      F95_check_sktrf_c(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrf_c(N, A, ipiv, "U", "P");
	F95_check_sktrf_c(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrf_c(N, A, ipiv, "L", "P");
	F95_check_sktrf_c(&N, Acpy, A, ipiv, "L", "P");
      }

      #ifdef __cplusplus

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrf(N, A, ipiv, "U", "N");
      F95_check_sktrf_c(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrf(N, A, ipiv, "L", "N");
      F95_check_sktrf_c(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrf(N, A, ipiv, "U", "P");
	F95_check_sktrf_c(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrf(N, A, ipiv, "L", "P");
	F95_check_sktrf_c(&N, Acpy, A, ipiv, "L", "P");
      }

      #endif

      /*check sktrd_x for various parameters*/
      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrd_c(N, A, tau, "U", "N");
      F95_check_sktrd_c(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrd_c(N, A, tau, "L", "N");
      F95_check_sktrd_c(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrd_c(N, A, tau, "U", "P");
	F95_check_sktrd_c(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrd_c(N, A, tau, "L", "P");
	F95_check_sktrd_c(&N, Acpy, A, tau, "L", "P");
      }

      #ifdef __cplusplus
      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrd(N, A, tau, "U", "N");
      F95_check_sktrd_c(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_c(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

      sktrd(N, A, tau, "L", "N");
      F95_check_sktrd_c(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrd(N, A, tau, "U", "P");
	F95_check_sktrd_c(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_c(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*N);

	sktrd(N, A, tau, "L", "P");
	F95_check_sktrd_c(&N, Acpy, A, tau, "L", "P");
      }

      #endif

      free(ipiv); free(tau); free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	floatcmplx *A, *Acpy;
	floatcmplx detQ, *Q;
	int i;

	A = (floatcmplx*)malloc(sizeof(floatcmplx)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (floatcmplx*)malloc(sizeof(floatcmplx)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}
	Q = (floatcmplx*)malloc(sizeof(floatcmplx)*N*N);
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(Acpy); free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd_c(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd_c(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd_c(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd_c(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	skbtrd_c(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	skbtrd_c(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd_c(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd_c(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd_c(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd_c(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	  skbtrd_c(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	  skbtrd_c(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#ifdef __cplusplus
	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_c(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_c(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(floatcmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = floatcmplx_one;

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_c(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#endif

	free(Q); free(Acpy); free(A);
      }
    }
  }
}

void check_decomp_z()
{
  int N, KD, idens;

  double densities[] = {0.1, 0.25, 0.5, 1.0};

  for(idens=0; idens<4; idens++) {
    /*dense matrices first*/

    for(N=1; N<=MAX_SIZE; N++) {
      doublecmplx *A, *Acpy;
      doublecmplx *tau;
      int *ipiv;

      A = (doublecmplx*)malloc(sizeof(doublecmplx)*N*N);
      if(!A) {
	printf("Failed to allocate memory!\n");
	exit(10);
      }
      Acpy = (doublecmplx*)malloc(sizeof(doublecmplx)*N*N);
      if(!Acpy) {
	printf("Failed to allocate memory!\n");
	free(A);
	exit(10);
      }
      tau = (doublecmplx*)malloc(sizeof(doublecmplx)*(N-1));
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(Acpy); free(A);
	exit(10);
      }
      ipiv = (int *)malloc(sizeof(int)*N);
      if(!tau) {
	printf("Failed to allocate memory!\n");
	free(tau), free(Acpy); free(A);
	exit(10);
      }

      /*check sktrf_x for various parameters*/
      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrf_z(N, A, ipiv, "U", "N");
      F95_check_sktrf_z(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrf_z(N, A, ipiv, "L", "N");
      F95_check_sktrf_z(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrf_z(N, A, ipiv, "U", "P");
	F95_check_sktrf_z(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrf_z(N, A, ipiv, "L", "P");
	F95_check_sktrf_z(&N, Acpy, A, ipiv, "L", "P");
      }

      #ifdef __cplusplus

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrf(N, A, ipiv, "U", "N");
      F95_check_sktrf_z(&N, Acpy, A, ipiv, "U", "N");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrf(N, A, ipiv, "L", "N");
      F95_check_sktrf_z(&N, Acpy, A, ipiv, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrf(N, A, ipiv, "U", "P");
	F95_check_sktrf_z(&N, Acpy, A, ipiv, "U", "P");

	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrf(N, A, ipiv, "L", "P");
	F95_check_sktrf_z(&N, Acpy, A, ipiv, "L", "P");
      }

      #endif

      /*check sktrd_x for various parameters*/
      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrd_z(N, A, tau, "U", "N");
      F95_check_sktrd_z(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrd_z(N, A, tau, "L", "N");
      F95_check_sktrd_z(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrd_z(N, A, tau, "U", "P");
	F95_check_sktrd_z(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrd_z(N, A, tau, "L", "P");
	F95_check_sktrd_z(&N, Acpy, A, tau, "L", "P");
      }

      #ifdef __cplusplus
      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrd(N, A, tau, "U", "N");
      F95_check_sktrd_z(&N, Acpy, A, tau, "U", "N");

      make_skew_mat_z(N, A, densities[idens]);
      memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

      sktrd(N, A, tau, "L", "N");
      F95_check_sktrd_z(&N, Acpy, A, tau, "L", "N");

      /*partial decompositions only for even N*/
      if(N%2 == 0) {
	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrd(N, A, tau, "U", "P");
	F95_check_sktrd_z(&N, Acpy, A, tau, "U", "P");

	make_skew_mat_z(N, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*N);

	sktrd(N, A, tau, "L", "P");
	F95_check_sktrd_z(&N, Acpy, A, tau, "L", "P");
      }

      #endif

      free(ipiv); free(tau); free(Acpy); free(A);
    }

    /*then banded matrices*/
    for(N=1; N<=MAX_SIZE; N++) {
      for(KD=0; KD<N; KD++) {
	doublecmplx *A, *Acpy;
	doublecmplx detQ, *Q;
	int i;

	A = (doublecmplx*)malloc(sizeof(doublecmplx)*N*(KD+1));
	if(!A) {
	  printf("Failed to allocate memory!\n");
	  exit(10);
	}
	Acpy = (doublecmplx*)malloc(sizeof(doublecmplx)*N*(KD+1));
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(A);
	  exit(10);
	}
	Q = (doublecmplx*)malloc(sizeof(doublecmplx)*N*N);
	if(!Acpy) {
	  printf("Failed to allocate memory!\n");
	  free(Acpy); free(A);
	  exit(10);
	}

	/*check skbpfa_x for various parameters*/
	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd_z(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd_z(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd_z(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd_z(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	skbtrd_z(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	skbtrd_z(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd_z(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd_z(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd_z(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd_z(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	  skbtrd_z(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	  skbtrd_z(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#ifdef __cplusplus

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "N");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "N");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "N");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "N");

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "V");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "V");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "V");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "V");

	make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	skbtrd(N, KD, A, &detQ, Q, "U", "N", "U");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "N", "U");

	make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	skbtrd(N, KD, A, &detQ, Q, "L", "N", "U");
	F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "N", "U");

	/*partial decompositions only for even N*/
	if(N%2 == 0) {
	  make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "N");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "N");

	  make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "N");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "N");

	  make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "V");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "V");

	  make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "V");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "V");

	  make_skew_mat_banded_z(N, KD, A, NULL, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	  skbtrd(N, KD, A, &detQ, Q, "U", "P", "U");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "U", "P", "U");

	  make_skew_mat_banded_z(N, KD, NULL, A, densities[idens]);
	  memcpy(Acpy, A, sizeof(doublecmplx)*N*(KD+1));
	  for(i=0; i<N*N; i++) Q[i] = doublecmplx_one;

	  skbtrd(N, KD, A, &detQ, Q, "L", "P", "U");
	  F95_check_skbtrd_z(&N, &KD, Acpy, A, &detQ, Q, "L", "P", "U");
	}

	#endif

	free(Q); free(Acpy); free(A);
      }
    }
  }
}


int main()
{
  printf("Checking interface for pfaffian calculation\n");

  check_pfaff_s();
  check_pfaff_d();
  check_pfaff_c();
  check_pfaff_z();

  printf("OK\n");

  printf("Checking interface for decompositions\n");

  check_decomp_s();
  check_decomp_d();
  check_decomp_c();
  check_decomp_z();

  printf("OK\n");

  printf("\nAll tests finished succesfully!\n");

  return 0;
}
