#ifndef PFAPACK_FORTRAN_H
#define PFAPACK_FORTRAN_H

#include "fortran.h"

/*****************************************************************
 Helper functions/macros and definitions to set up dense or banded
 matrices in Fortran format
******************************************************************/

#ifdef __cplusplus

#include <algorithm>
#include <cstdio>
#include <cstdlib>

inline int dense_fortran(int i, int j, int ldim)
{
  return i-1+(j-1)*ldim;
}

inline int bandlower_fortran(int i, int j, int N, int KD)
{
  if(j<=i && i<=std::min(N,j+KD)) {
    return i-j + (KD+1)*(j-1);
  }
  else {
    std::printf("Index out of range for lower band matrix!\n");
    std::exit(10);
  }
}

inline int bandupper_fortran(int i, int j, int N, int KD)
{
  if(i<=j && i>=std::max(1,j-KD)) {
    return KD+i-j + (KD+1)*(j-1);
  }
  else {
    std::printf("Index out of range for upper band matrix!\n");
    std::exit(10);
  }
}

#else

#if __STDC_VERSION__ >= 199901L

#include <stdio.h>
#include <stdlib.h>

inline int dense_fortran(int i, int j, int ldim)
{
  return i-1+(j-1)*ldim;
}

inline int bandlower_fortran(int i, int j, int N, int KD)
{
  if(j<=i && i<=(N<j+KD ? N : j+KD)) {
    return i-j + (KD+1)*(j-1);
  }
  else {
    printf("Index out of range for lower band matrix!");
    exit(10);
  }
}

inline int bandupper_fortran(int i, int j, int N, int KD)
{
  if(i<=j && i>=(1>j-KD ? 1 : j-KD)) {
    return KD+i-j + (KD+1)*(j-1);
  }
  else {
    printf("Index out of range for upper band matrix!");
    exit(10);
  }
}

#else

#define dense_fortran(i, j, ldim) (i-1+(j-1)*ldim)
#define bandlower_fortran(i, j, N, KD) (i-j + (KD+1)*(j-1))
#define bandupper_fortran(i, j, N, KD) (KD+i-j + (KD+1)*(j-1))

#endif

#endif

#define PFAPACK_ssktrf fortran_name(ssktrf, SSKTRF)
#define PFAPACK_dsktrf fortran_name(dsktrf, DSKTRF)
#define PFAPACK_csktrf fortran_name(csktrf, CSKTRF)
#define PFAPACK_zsktrf fortran_name(zsktrf, ZSKTRF)

#define PFAPACK_ssktrd fortran_name(ssktrd, SSKTRD)
#define PFAPACK_dsktrd fortran_name(dsktrd, DSKTRD)
#define PFAPACK_csktrd fortran_name(csktrd, CSKTRD)
#define PFAPACK_zsktrd fortran_name(zsktrd, ZSKTRD)

#define PFAPACK_sskbtrd fortran_name(sskbtrd, SSKTRD)
#define PFAPACK_dskbtrd fortran_name(dskbtrd, DSKTRD)
#define PFAPACK_cskbtrd fortran_name(cskbtrd, CSKTRD)
#define PFAPACK_zskbtrd fortran_name(zskbtrd, ZSKTRD)

#define PFAPACK_sskpfa fortran_name(sskpfa, SSKPFA)
#define PFAPACK_dskpfa fortran_name(dskpfa, DSKPFA)
#define PFAPACK_cskpfa fortran_name(cskpfa, CSKPFA)
#define PFAPACK_zskpfa fortran_name(zskpfa, ZSKPFA)

#define PFAPACK_sskpf10 fortran_name(sskpf10, SSKPF10)
#define PFAPACK_dskpf10 fortran_name(dskpf10, DSKPF10)
#define PFAPACK_cskpf10 fortran_name(cskpf10, CSKPF10)
#define PFAPACK_zskpf10 fortran_name(zskpf10, ZSKPF10)

#define PFAPACK_sskbpfa fortran_name(sskbpfa, SSKBPFA)
#define PFAPACK_dskbpfa fortran_name(dskbpfa, DSKBPFA)
#define PFAPACK_cskbpfa fortran_name(cskbpfa, CSKBPFA)
#define PFAPACK_zskbpfa fortran_name(zskbpfa, ZSKBPFA)

#define PFAPACK_sskbpf10 fortran_name(sskbpf10, SSKBPF10)
#define PFAPACK_dskbpf10 fortran_name(dskbpf10, DSKBPF10)
#define PFAPACK_cskbpf10 fortran_name(cskbpf10, CSKBPF10)
#define PFAPACK_zskbpf10 fortran_name(zskbpf10, ZSKBPF10)


#ifdef __cplusplus
extern "C" {
#endif

/* prototypes for the PFAPACK routines S,D,C,ZSKTRF */

void PFAPACK_ssktrf(const char *, const char *, const int *,
		    float *, const int *, int *, float *,
		    const int *, int *);

void PFAPACK_dsktrf(const char *, const char *, const int *,
		    double *, const int *, int *, double *,
		    const int *, int *);

void PFAPACK_csktrf(const char *, const char *, const int *,
		    floatcmplx *, const int *, int *,
		    floatcmplx *, const int *, int *);

void PFAPACK_zsktrf(const char *, const char *, const int *,
		    doublecmplx *, const int *, int *,
		    doublecmplx *, const int *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKTF2 */

void PFAPACK_ssktf2(const char *, const char *, const int *,
		    float *, const int *, int *, int *);

void PFAPACK_dsktf2(const char *, const char *, const int *,
		    double *, const int *, int *, int *);

void PFAPACK_csktf2(const char *, const char *, const int *,
		    floatcmplx *, const int *, int *,
		    int *);

void PFAPACK_zsktf2(const char *, const char *, const int *,
		    doublecmplx *, const int *, int *,
		    int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKTRD */

void PFAPACK_ssktrd(const char *, const char *, const int *,
		    float *, const int *, float *, float *,
		    float *, const int *, int *);

void PFAPACK_dsktrd(const char *, const char *, const int *,
		    double *, const int *, double *, double *,
		    double *, const int *, int *);

void PFAPACK_csktrd(const char *, const char *, const int *,
		    floatcmplx *, const int *, float *,
		    floatcmplx *,
		    floatcmplx *, const int *, int *);

void PFAPACK_zsktrd(const char *, const char *, const int *,
		    doublecmplx *, const int *,
		    double *, doublecmplx *,
		    doublecmplx *, const int *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKTD2 */

void PFAPACK_ssktd2(const char *, const char *, const int *,
		    float *, const int *, float *, float *, int *);

void PFAPACK_dsktd2(const char *, const char *, const int *,
		    double *, const int *, double *, double *,
		    int *);

void PFAPACK_csktd2(const char *, const char *, const int *,
		    floatcmplx *, const int *, float *,
		    floatcmplx *, int *);

void PFAPACK_zsktd2(const char *, const char *, const int *,
		    doublecmplx *, const int *,
		    double *, doublecmplx *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKPFA */

void PFAPACK_sskpfa(const char *, const char *, const int *,
		    float *, const int *, float *, int *,
		    float *, const int *, int *);

void PFAPACK_dskpfa(const char *, const char *, const int *,
		    double *, const int *, double *, int *,
		    double *, const int *, int *);

void PFAPACK_cskpfa(const char *, const char *, const int *,
		    floatcmplx *, const int *,
		    floatcmplx *, int *,
		    floatcmplx *, const int *, float *,
		    int *);

void PFAPACK_zskpfa(const char *, const char *, const int *,
		    doublecmplx *, const int *,
		    doublecmplx *, int *,
		    doublecmplx *, const int *,
		    double *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKPF10 */

void PFAPACK_sskpf10(const char *, const char *, const int *,
		     float *, const int *, float *, int *,
		     float *, const int *, int *);

void PFAPACK_dskpf10(const char *, const char *, const int *,
		     double *, const int *, double *, int *,
		     double *, const int *, int *);

void PFAPACK_cskpf10(const char *, const char *, const int *,
		     floatcmplx *, const int *,
		     floatcmplx *, int *,
		     floatcmplx *, const int *,
		     float *, int *);

void PFAPACK_zskpf10(const char *, const char *, const int *,
		     doublecmplx *, const int *,
		     doublecmplx *, int *,
		     doublecmplx *, const int *,
		     double *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKBTRD */

void PFAPACK_sskbtrd(const char *, const char *, const char *,
		     const int *, const int *, float *,
		     const int *, float *, float *, const int *,
		     float *, int *);

void PFAPACK_dskbtrd(const char *, const char *, const char *,
		     const int *, const int *, double *,
		     const int *, double *, double *,
		     const int *, double *, int *);

void PFAPACK_cskbtrd(const char *, const char *, const char *,
		     const int *, const int *,
		     floatcmplx *, const int *,
		     float *, floatcmplx *,
		     floatcmplx *, const int *,
		     floatcmplx *, float *, int *);

void PFAPACK_zskbtrd(const char *, const char *, const char *,
		     const int *, const int *,
		     doublecmplx *, const int *,
		     double *, doublecmplx *,
		     doublecmplx *, const int *,
		     doublecmplx *, double *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKBPFA */

void PFAPACK_sskbpfa(const char *, const int *, const int *,
		     float *, const int *, float *, float *,
		     int *);

void PFAPACK_dskbpfa(const char *, const int *, const int *,
		     double *, const int *, double *, double *,
		     int *);

void PFAPACK_cskbpfa(const char *, const int *, const int *,
		     floatcmplx *, const int *,
		     floatcmplx *,
		     floatcmplx *, float *, int *);

void PFAPACK_zskbpfa(const char *, const int *, const int *,
		     doublecmplx *, const int *,
		     doublecmplx *,
		     doublecmplx *, double *, int *);

/* prototypes for the PFAPACK routines S,D,C,ZSKBPF10 */

void PFAPACK_sskbpf10(const char *, const int *, const int *,
		      float *, const int *, float *, float *,
		      int *);

void PFAPACK_dskbpf10(const char *, const int *, const int *,
		      double *, const int *, double *, double *,
		      int *);

void PFAPACK_cskbpf10(const char *, const int *, const int *,
		      floatcmplx *, const int *,
		      floatcmplx *,
		      floatcmplx *, float *, int *);

void PFAPACK_zskbpf10(const char *, const int *, const int *,
		      doublecmplx *, const int *,
		      doublecmplx *,
		      doublecmplx *, double *, int *);

#ifdef __cplusplus
}
#endif

#endif
