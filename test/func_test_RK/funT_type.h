/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: Macro define to realize 4 data type.
 * Author: linhouzhong
 * Create: 2022-06-26
 */

#pragma once

#include <complex.h>

#ifdef SINGLE
#define TestTypeName SINGLE
typedef float funcTdataAccu;
typedef float funcTdataType;
#define _sysv_rk_ ssysv_rk_
#define _CreatMatrix sCreatMatrix
#define _FuncTest sFuncTest
#define _pot02_ spot02_
#define _get04_ sget04_
#define _gemm_ sgemm_
#define _lansy_ slansy_
#define _sycon_ ssycon_
#define _sytrf_ ssytrf_
#endif
#ifdef DOUBLE
#define TestTypeName DOUBLE
typedef double funcTdataAccu;
typedef double funcTdataType;
#define _sysv_rk_ dsysv_rk_
#define _CreatMatrix dCreatMatrix
#define _FuncTest dFuncTest
#define _pot02_ dpot02_
#define _get04_ dget04_
#define _gemm_ dgemm_
#define _lansy_ dlansy_
#define _sycon_ dsycon_
#define _sytrf_ dsytrf_
#endif
#ifdef COMPLEX
#define TestTypeName COMPLEX
typedef float funcTdataAccu;
typedef float complex funcTdataType;
#define _sysv_rk_ csysv_rk_
#define _CreatMatrix cCreatMatrix
#define _FuncTest cFuncTest
#define _pot02_ csyt02_
#define _get04_ cget04_
#define _gemm_ cgemm_
#define _lansy_ clansy_
#define _sycon_ csycon_
#define _sytrf_ csytrf_
#endif
#ifdef COMPLEX16
#define TestTypeName COMPLEX16
typedef double funcTdataAccu;
typedef double complex funcTdataType;
#define _sysv_rk_ zsysv_rk_
#define _CreatMatrix zCreatMatrix
#define _FuncTest zFuncTest
#define _pot02_ zsyt02_
#define _get04_ zget04_
#define _gemm_ zgemm_
#define _lansy_ zlansy_
#define _sycon_ zsycon_
#define _sytrf_ zsytrf_
#endif

