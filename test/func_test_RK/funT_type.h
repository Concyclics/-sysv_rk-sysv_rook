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
#endif
#ifdef DOUBLE
#define TestTypeName DOUBLE
typedef double funcTdataAccu;
typedef double funcTdataType;
#define _sysv_rk_ dsysv_rk_
#define _CreatMatrix dCreatMatrix
#define _FuncTest dFuncTest
#define _pot02_ dpot02_
#endif
#ifdef COMPLEX
#define TestTypeName COMPLEX
typedef float funcTdataAccu;
typedef float complex funcTdataType;
#define _sysv_rk_ csysv_rk_
#define _CreatMatrix cCreatMatrix
#define _FuncTest cFuncTest
#define _pot02_ csyt02_
#endif
#ifdef COMPLEX16
#define TestTypeName COMPLEX16
typedef double funcTdataAccu;
typedef double complex funcTdataType;
#define _sysv_rk_ zsysv_rk_
#define _CreatMatrix zCreatMatrix
#define _FuncTest zFuncTest
#define _pot02_ zsyt02_
#endif
