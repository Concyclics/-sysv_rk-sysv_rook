#pragma once

#include <complex.h>

#ifdef SINGLE
typedef float funcTdataAccu;
typedef float funcTdataType;
#define _sysv_rook_ ssysv_rook_
#define _CreatMatrix sCreatMatrixMy
#define _FuncTest sFuncTest
#define _pot02_ spot02_
#endif
#ifdef DOUBLE
typedef double funcTdataAccu;
typedef double funcTdataType;
#define _sysv_rook_ dsysv_rook_
#define _CreatMatrix dCreatMatrixMy
#define _FuncTest dFuncTest
#define _pot02_ dpot02_
#endif
#ifdef COMPLEX
typedef float funcTdataAccu;
typedef float complex funcTdataType;
#define _sysv_rook_ csysv_rook_
#define _CreatMatrix cCreatMatrixMy
#define _FuncTest cFuncTest
#define _pot02_ cpot02_
#endif
#ifdef COMPLEX16
typedef double funcTdataAccu;
typedef double complex funcTdataType;
#define _sysv_rook_ zsysv_rook_
#define _CreatMatrix zCreatMatrixMy
#define _FuncTest zFuncTest
#define _pot02_ zpot02_
#endif
