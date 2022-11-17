/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: Macro define to realize 4 data type.
 * Author: linhouzhong
 * Create: 2022-06-26
 */

#pragma once

#include <complex.h>

#ifdef SINGLE
typedef float funcTdataAccu;
typedef float funcTdataType;
#define _sysv_rook_ ssysv_rook_
#define _CreatMatrix sCreatMatrix
#define _FuncTest sFuncTest
#define _pot02_ spot02_
#define _get04_ sget04_
#define _gemm_ sgemm_
#endif
#ifdef DOUBLE
typedef double funcTdataAccu;
typedef double funcTdataType;
#define _sysv_rook_ dsysv_rook_
#define _CreatMatrix dCreatMatrix
#define _FuncTest dFuncTest
#define _pot02_ dpot02_
#define _get04_ dget04_
#define _gemm_ dgemm_
#endif
#ifdef COMPLEX
typedef float funcTdataAccu;
typedef float complex funcTdataType;
#define _sysv_rook_ csysv_rook_
#define _CreatMatrix cCreatMatrix
#define _FuncTest cFuncTest
#define _pot02_ csyt02_
#define _get04_ cget04_
#define _gemm_ cgemm_
#endif
#ifdef COMPLEX16
typedef double funcTdataAccu;
typedef double complex funcTdataType;
#define _sysv_rook_ zsysv_rook_
#define _CreatMatrix zCreatMatrix
#define _FuncTest zFuncTest
#define _pot02_ zsyt02_
#define _get04_ zget04_
#define _gemm_ zgemm_
#endif
