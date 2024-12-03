#include "fortran.h"

#define F95_check_sktrf_s fortran_name(check_sktrf_s, CHECK_SKTRF_S)
#define F95_check_sktrf_d fortran_name(check_sktrf_d, CHECK_SKTRF_D)
#define F95_check_sktrf_c fortran_name(check_sktrf_c, CHECK_SKTRF_C)
#define F95_check_sktrf_z fortran_name(check_sktrf_z, CHECK_SKTRF_Z)

#define F95_check_sktrd_s fortran_name(check_sktrd_s, CHECK_SKTRD_S)
#define F95_check_sktrd_d fortran_name(check_sktrd_d, CHECK_SKTRD_D)
#define F95_check_sktrd_c fortran_name(check_sktrd_c, CHECK_SKTRD_C)
#define F95_check_sktrd_z fortran_name(check_sktrd_z, CHECK_SKTRD_Z)

#define F95_check_skbtrd_s fortran_name(check_skbtrd_s, CHECK_SKBTRD_S)
#define F95_check_skbtrd_d fortran_name(check_skbtrd_d, CHECK_SKBTRD_D)
#define F95_check_skbtrd_c fortran_name(check_skbtrd_c, CHECK_SKBTRD_C)
#define F95_check_skbtrd_z fortran_name(check_skbtrd_z, CHECK_SKBTRD_Z)

#define F95_check_skpfa_s fortran_name(check_skpfa_s, CHECK_SKPFA_S)
#define F95_check_skpfa_d fortran_name(check_skpfa_d, CHECK_SKPFA_D)
#define F95_check_skpfa_c fortran_name(check_skpfa_c, CHECK_SKPFA_C)
#define F95_check_skpfa_z fortran_name(check_skpfa_z, CHECK_SKPFA_Z)

#define F95_check_skpf10_s fortran_name(check_skpf10_s, CHECK_SKPF10_S)
#define F95_check_skpf10_d fortran_name(check_skpf10_d, CHECK_SKPF10_D)
#define F95_check_skpf10_c fortran_name(check_skpf10_c, CHECK_SKPF10_C)
#define F95_check_skpf10_z fortran_name(check_skpf10_z, CHECK_SKPF10_Z)

#define F95_check_skbpfa_s fortran_name(check_skbpfa_s, CHECK_SKBPFA_S)
#define F95_check_skbpfa_d fortran_name(check_skbpfa_d, CHECK_SKBPFA_D)
#define F95_check_skbpfa_c fortran_name(check_skbpfa_c, CHECK_SKBPFA_C)
#define F95_check_skbpfa_z fortran_name(check_skbpfa_z, CHECK_SKBPFA_Z)

#define F95_check_skbpf10_s fortran_name(check_skbpf10_s, CHECK_SKBPF10_S)
#define F95_check_skbpf10_d fortran_name(check_skbpf10_d, CHECK_SKBPF10_D)
#define F95_check_skbpf10_c fortran_name(check_skbpf10_c, CHECK_SKBPF10_C)
#define F95_check_skbpf10_z fortran_name(check_skbpf10_z, CHECK_SKBPF10_Z)

#ifdef __cplusplus
extern "C" {
#endif

void F95_check_skpfa_s(int *, float *, float *, float *, const char *,
		       const char *);

void F95_check_skpf10_s(int *, float *, float *, float *, const char *,
			const char *);

void F95_check_skbpfa_s(int *, int *, float *, float *, float *,
			const char *);

void F95_check_skbpf10_s(int *, int *, float *, float *, float *,
			 const char *);

void F95_check_skpfa_d(int *, double *, double *, double *, const char *,
		       const char *);

void F95_check_skpf10_d(int *, double *, double *, double *, const char *,
			const char *);

void F95_check_skbpfa_d(int *, int *, double *, double *, double *,
			const char *);

void F95_check_skbpf10_d(int *, int *, double *, double *, double *,
			 const char *);

void F95_check_skpfa_c(int *, floatcmplx *, floatcmplx *, floatcmplx *,
		       const char *, const char *);

void F95_check_skpf10_c(int *, floatcmplx *, floatcmplx *, floatcmplx *,
			const char *, const char *);

void F95_check_skbpfa_c(int *, int *, floatcmplx *, floatcmplx *,
			floatcmplx *, const char *);

void F95_check_skbpf10_c(int *, int *, floatcmplx *, floatcmplx *,
			 floatcmplx *, const char *);

void F95_check_skpfa_z(int *, doublecmplx *, doublecmplx *, doublecmplx *,
		       const char *, const char *);

void F95_check_skpf10_z(int *, doublecmplx *, doublecmplx *,
			doublecmplx *, const char *, const char *);

void F95_check_skbpfa_z(int *, int *, doublecmplx *, doublecmplx *,
			doublecmplx *, const char *);

void F95_check_skbpf10_z(int *, int *, doublecmplx *, doublecmplx *,
			 doublecmplx *, const char *);

void F95_check_sktrf_s(int *, float *, float *, int *, const char *,
		       const char *);

void F95_check_sktrd_s(int *, float *, float *, float *, const char *,
		       const char *);

void F95_check_skbtrd_s(int *, int *, float *, float *, float *, float *,
			const char *, const char *, const char *);

void F95_check_sktrf_d(int *, double *, double *, int *, const char *,
		       const char *);

void F95_check_sktrd_d(int *, double *, double *, double *, const char *,
		       const char *);

void F95_check_skbtrd_d(int *, int *, double *, double *, double *,
			double *, const char *, const char *,
			const char *);

void F95_check_sktrf_c(int *, floatcmplx *, floatcmplx *, int *,
		       const char *, const char *);

void F95_check_sktrd_c(int *, floatcmplx *, floatcmplx *, floatcmplx *,
		       const char *, const char *);

void F95_check_skbtrd_c(int *, int *, floatcmplx *, floatcmplx *,
			floatcmplx *, floatcmplx *, const char *,
			const char *, const char *);

void F95_check_sktrf_z(int *, doublecmplx *, doublecmplx *, int *,
		       const char *, const char *);

void F95_check_sktrd_z(int *, doublecmplx *, doublecmplx *, doublecmplx *,
		       const char *, const char *);

void F95_check_skbtrd_z(int *, int *, doublecmplx *, doublecmplx *,
			doublecmplx *, doublecmplx *, const char *,
			const char *, const char *);

#ifdef __cplusplus
}
#endif
