#ifndef PFAPACK_H
#define PFAPACK_H

#include "fortran_pfapack.h"

#ifdef __cplusplus
extern "C" {
#endif

int skpfa_s(int, float *, float *, const char *, const char *);
int skpfa_d(int, double *, double *, const char *, const char *);
int skpfa_c(int, floatcmplx *, floatcmplx *,
	    const char *, const char *);
int skpfa_z(int, doublecmplx *, doublecmplx *,
	    const char *, const char *);

int skpf10_s(int, float *, float *, const char *, const char *);
int skpf10_d(int, double *, double *, const char *, const char *);
int skpf10_c(int, floatcmplx *, floatcmplx *,
	     const char *, const char *);
int skpf10_z(int, doublecmplx *, doublecmplx *,
	     const char *, const char *);

int skbpfa_s(int, int, float *, float *, const char *);
int skbpfa_d(int, int, double *, double *, const char *);
int skbpfa_c(int, int, floatcmplx *, floatcmplx *,
	     const char *);
int skbpfa_z(int, int, doublecmplx *, doublecmplx *,
	     const char *);

int skbpf10_s(int, int, float *, float *, const char *);
int skbpf10_d(int, int, double *, double *, const char *);
int skbpf10_c(int, int, floatcmplx *, floatcmplx *,
	      const char *);
int skbpf10_z(int, int, doublecmplx *, doublecmplx *,
	      const char *);

int sktrf_s(int, float *, int *, const char *, const char *);
int sktrf_d(int, double *, int *, const char *, const char *);
int sktrf_c(int, floatcmplx *, int *, const char *, const char *);
int sktrf_z(int, doublecmplx *, int *, const char *, const char *);

int sktrd_s(int, float *, float *, const char *, const char *);
int sktrd_d(int, double *, double *, const char *, const char *);
int sktrd_c(int, floatcmplx *, floatcmplx *, const char *, const char *);
int sktrd_z(int, doublecmplx *, doublecmplx *, const char *, const char *);

int skbtrd_s(int, int, float *, float *, float *,
	     const char *, const char *, const char *);
int skbtrd_d(int, int, double*, double *, double *,
	     const char *, const char *, const char *);
int skbtrd_c(int, int, floatcmplx *, floatcmplx *, floatcmplx *,
	     const char *, const char *, const char *);
int skbtrd_z(int, int, doublecmplx *, doublecmplx *, doublecmplx *,
	     const char *, const char *, const char *);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* overloaded versions of skpfa_x */
inline int skpfa(int N, float *A, float *PFAFF,
		 const char *UPLO, const char *MTHD)
{ return skpfa_s(N, A, PFAFF, UPLO, MTHD); }

inline int skpfa(int N, double *A, double *PFAFF,
		 const char *UPLO, const char *MTHD)
{ return skpfa_d(N, A, PFAFF, UPLO, MTHD); }

inline int skpfa(int N, floatcmplx *A, floatcmplx *PFAFF,
		 const char *UPLO, const char *MTHD)
{ return skpfa_c(N, A, PFAFF, UPLO, MTHD); }

inline int skpfa(int N, doublecmplx *A, doublecmplx *PFAFF,
		 const char *UPLO, const char *MTHD)
{ return skpfa_z(N, A, PFAFF, UPLO, MTHD); }

/*overloaded versions of skpf10_x */
inline int skpf10(int N, float *A, float *PFAFF,
		  const char *UPLO, const char *MTHD)
{ return skpf10_s(N, A, PFAFF, UPLO, MTHD); }

inline int skpf10(int N, double *A, double *PFAFF,
		  const char *UPLO, const char *MTHD)
{ return skpf10_d(N, A, PFAFF, UPLO, MTHD); }

inline int skpf10(int N, floatcmplx *A, floatcmplx *PFAFF,
		  const char *UPLO, const char *MTHD)
{ return skpf10_c(N, A, PFAFF, UPLO, MTHD); }

inline int skpf10(int N, doublecmplx *A, doublecmplx *PFAFF,
		  const char *UPLO, const char *MTHD)
{ return skpf10_z(N, A, PFAFF, UPLO, MTHD); }

/* overloaded versions of skbpfa_x */
inline int skbpfa(int N, int KD, float *A, float *PFAFF,
		  const char *UPLO)
{ return skbpfa_s(N, KD, A, PFAFF, UPLO); }

inline int skbpfa(int N, int KD, double *A, double *PFAFF,
		  const char *UPLO)
{ return skbpfa_d(N, KD, A, PFAFF, UPLO); }

inline int skbpfa(int N, int KD, floatcmplx *A, floatcmplx *PFAFF,
		  const char *UPLO)
{ return skbpfa_c(N, KD, A, PFAFF, UPLO); }

inline int skbpfa(int N, int KD, doublecmplx *A, doublecmplx *PFAFF,
		  const char *UPLO)
{ return skbpfa_z(N, KD, A, PFAFF, UPLO); }

/* overloaded versions of skbpf10_x */
inline int skbpf10(int N, int KD, float *A, float *PFAFF,
		   const char *UPLO)
{ return skbpf10_s(N, KD, A, PFAFF, UPLO); }

inline int skbpf10(int N, int KD, double *A, double *PFAFF,
		   const char *UPLO)
{ return skbpf10_d(N, KD, A, PFAFF, UPLO); }

inline int skbpf10(int N, int KD, floatcmplx *A, floatcmplx *PFAFF,
		   const char *UPLO)
{ return skbpf10_c(N, KD, A, PFAFF, UPLO); }

inline int skbpf10(int N, int KD, doublecmplx *A, doublecmplx *PFAFF,
		   const char *UPLO)
{ return skbpf10_z(N, KD, A, PFAFF, UPLO); }

/* overloaded versions of sktrf_x */
inline int sktrf(int N, float *A, int *IPIV, const char *UPLO,
		 const char *MODE)
{ return  sktrf_s(N, A, IPIV, UPLO, MODE); }

inline int sktrf(int N, double *A, int *IPIV, const char *UPLO,
		 const char *MODE)
{ return  sktrf_d(N, A, IPIV, UPLO, MODE); }

inline int sktrf(int N, floatcmplx *A, int *IPIV, const char *UPLO,
		 const char *MODE)
{ return  sktrf_c(N, A, IPIV, UPLO, MODE); }

inline int sktrf(int N, doublecmplx *A, int *IPIV, const char *UPLO,
		 const char *MODE)
{ return  sktrf_z(N, A, IPIV, UPLO, MODE); }

/* overloaded versions of sktrd_x */
inline int sktrd(int N, float *A, float *TAU, const char *UPLO,
		 const char *MODE)
{ return sktrd_s(N, A, TAU, UPLO, MODE); }

inline int sktrd(int N, double *A, double *TAU, const char *UPLO,
		 const char *MODE)
{ return sktrd_d(N, A, TAU, UPLO, MODE); }

inline int sktrd(int N, floatcmplx *A, floatcmplx *TAU, const char *UPLO,
		 const char *MODE)
{ return sktrd_c(N, A, TAU, UPLO, MODE); }

inline int sktrd(int N, doublecmplx *A, doublecmplx *TAU, const char *UPLO,
		 const char *MODE)
{ return sktrd_z(N, A, TAU, UPLO, MODE); }

/* overloaded versions of skbtrd_x */
inline int skbtrd(int N, int KD, float *A, float *DETQ, float *Q,
		  const char *UPLO, const char *MODE, const char *VECT)
{ return skbtrd_s(N, KD, A, DETQ, Q, UPLO, MODE, VECT); }

inline int skbtrd(int N, int KD, double *A, double *DETQ, double *Q,
		  const char *UPLO, const char *MODE, const char *VECT)
{ return skbtrd_d(N, KD, A, DETQ, Q, UPLO, MODE, VECT); }

inline int skbtrd(int N, int KD, floatcmplx *A, floatcmplx *DETQ,
		  floatcmplx *Q, const char *UPLO, const char *MODE,
		  const char *VECT)
{ return skbtrd_c(N, KD, A, DETQ, Q, UPLO, MODE, VECT); }

inline int skbtrd(int N, int KD, doublecmplx *A, doublecmplx *DETQ,
		  doublecmplx *Q, const char *UPLO, const char *MODE,
		  const char *VECT)
{ return skbtrd_z(N, KD, A, DETQ, Q, UPLO, MODE, VECT); }

#endif

#endif
