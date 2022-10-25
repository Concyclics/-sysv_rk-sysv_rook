#pragma once

#include <complex.h>

#pragma once

#ifdef SINGLE
#define SYTRF_NAME "SSYTRF"
#define SYTRS_NAME "SSYTRS"
#define _opla_ dopla_
#define _sysv_rk_ ssysv_rk_
#endif
#ifdef DOUBLE
#define SYTRF_NAME "DSYTRF"
#define SYTRS_NAME "DSYTRS"
#define _opla_ dopla_
#define _sysv_rk_ dsysv_rk_
#endif
#ifdef COMPLEX
#include <complex.h>
#define SYTRF_NAME "CSYTRF"
#define SYTRS_NAME "CSYTRS"
#define _opla_ dopla_
#define _sysv_rk_ csysv_rk_
#endif
#ifdef COMPLEX16
#include <complex.h>
#define SYTRF_NAME "ZSYTRF"
#define SYTRS_NAME "ZSYTRS"
#define _opla_ dopla_
#define _sysv_rk_ zsysv_rk_
#endif
