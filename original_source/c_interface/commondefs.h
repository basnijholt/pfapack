#ifdef __cplusplus

#include <cstdlib>
#include <cctype>
#include <cstdio>

using std::malloc;
using std::free;
using std::toupper;
using std::printf;
using std::exit;

#else

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#endif

#include "fortran_pfapack.h"

#ifdef CPLUSPLUS_COMPLEX

static const floatcmplx floatcmplx_one = std::complex<float>(1.0);
static const doublecmplx doublecmplx_one = std::complex<double>(1.0);

static const floatcmplx floatcmplx_zero = std::complex<float>(0.0);
static const doublecmplx doublecmplx_zero = std::complex<double>(0.0);

#elif defined(C99_COMPLEX)

static const floatcmplx floatcmplx_one = 1;
static const doublecmplx doublecmplx_one = 1;

static const floatcmplx floatcmplx_zero = 0;
static const doublecmplx doublecmplx_zero = 0;

#elif defined(STRUCT_COMPLEX)

static const floatcmplx floatcmplx_one = { 1, 0 };
static const doublecmplx doublecmplx_one = { 1, 0 };

static const floatcmplx floatcmplx_zero = { 0, 0 };
static const doublecmplx doublecmplx_zero = { 0, 0 };

#endif
