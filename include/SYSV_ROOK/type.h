/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: Macro define to realize 4 data type.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

#include <kml_service.h>
#include <omp.h>

#define MAX(x, y) KmlMax(x, y)
#define MIN(x, y) KmlMin(x, y)

#ifdef SINGLE
typedef float dataAccu;
typedef float dataType;
#define ZERO 0
#define ONE 1
#define EIGHT 8
#define SEVTEN 17
#define T_ONE ONE
#define P_ONE ONE
#define P_NEG_ONE -ONE
#define P_R1 r1
#define ABS_ fabs
#define lamch_ Slamch
#define SYTRF_ROOK ssytrf_rook_
#define SYTRS_ROOK ssytrs_rook_
#define SYTF2_ROOK ssytf2_rook_
#define LASYF_ROOK slasyf_rook_
#define SYSV_ROOK ssysv_rook_
#endif
#ifdef DOUBLE
typedef double dataAccu;
typedef double dataType;
#define ZERO 0
#define ONE 1
#define EIGHT 8
#define SEVTEN 17
#define T_ONE ONE
#define P_ONE ONE
#define P_NEG_ONE -ONE
#define P_R1 r1
#define ABS_ fabs
#define lamch_ Dlamch
#define SYTRF_ROOK dsytrf_rook_
#define SYTRS_ROOK dsytrs_rook_
#define SYTF2_ROOK dsytf2_rook_
#define LASYF_ROOK dlasyf_rook_
#define SYSV_ROOK dsysv_rook_
#endif
#ifdef COMPLEX
#include <complex.h>
typedef float dataAccu;
typedef kml_complex_float dataType;
#define ZERO 0
#define ONE 1
#define EIGHT 8
#define SEVTEN 17
#define T_ONE CONE
#define P_ONE &CONE
#define P_NEG_ONE &NEG_CONE
#define P_R1 &r1
#define ABS_ CABS
#define CABS(Z) (fabs(creal(Z)) + fabs(cimag(Z)))
#define lamch_ Slamch
#define SYTRF_ROOK csytrf_rook_
#define SYTRS_ROOK csytrs_rook_
#define SYTF2_ROOK csytf2_rook_
#define LASYF_ROOK clasyf_rook_
#define SYSV_ROOK csysv_rook_
#endif
#ifdef COMPLEX16
#include <complex.h>
typedef double dataAccu;
typedef kml_complex_double dataType;
#define ZERO 0
#define ONE 1
#define EIGHT 8
#define SEVTEN 17
#define T_ONE CONE
#define P_ONE &CONE
#define P_NEG_ONE &NEG_CONE
#define P_R1 &r1
#define ABS_ CABS
#define CABS(Z) (fabs(creal(Z)) + fabs(cimag(Z)))
#define lamch_ Dlamch
#define SYTRF_ROOK zsytrf_rook_
#define SYTRS_ROOK zsytrs_rook_
#define SYTF2_ROOK zsytf2_rook_
#define LASYF_ROOK zlasyf_rook_
#define SYSV_ROOK zsysv_rook_
#endif

#ifdef SINGLE
typedef float blasDataAccu;
typedef float blasDataType;
#define blasParamType float
#define I_AMAX isamax_
#define COPY_ scopy_
#define SCAL_ sscal_
#define SWAP_ sswap_
#define SYR_ ssyr_
#define GER_ sger_
#define GEMV_ sgemv_
#define GEMM_ sgemm_
#define TRSM_ strsm_
#endif

#ifdef DOUBLE
typedef double blasDataAccu;
typedef double blasDataType;
#define blasParamType double
#define I_AMAX idamax_
#define COPY_ dcopy_
#define SCAL_ dscal_
#define SWAP_ dswap_
#define SYR_ dsyr_
#define GER_ dger_
#define GEMV_ dgemv_
#define GEMM_ dgemm_
#define TRSM_ dtrsm_
#endif

#ifdef COMPLEX
typedef float blasDataAccu;
typedef kml_complex_float blasDataType;
#define blasParamType void*
#define I_AMAX icamax_
#define COPY_ ccopy_
#define SCAL_ cscal_
#define SWAP_ cswap_
#define SYR_ csyr_
#define GER_ cgeru_
#define GEMV_ cgemv_
#define GEMM_ cgemm_
#define TRSM_ ctrsm_
#endif

#ifdef COMPLEX16
typedef double blasDataAccu;
typedef kml_complex_double blasDataType;
#define blasParamType void*
#define I_AMAX izamax_
#define COPY_ zcopy_
#define SCAL_ zscal_
#define SWAP_ zswap_
#define SYR_ zsyr_
#define GER_ zgeru_
#define GEMV_ zgemv_
#define GEMM_ zgemm_
#define TRSM_ ztrsm_
#endif

// I_AMAX
int I_AMAX(const int* N, blasDataType* x, const int* incX);

// COPY
void COPY_(const int* n,
           const blasDataType* x,
           const int* incX,
           blasDataType* y,
           const int* incY);

// SCAL
void SCAL_(const int* N,
           const blasDataType* alpha,
           blasDataType* X,
           const int* incX);

// SWAP
void SWAP_(const int* n,
           blasDataType* x,
           const int* incX,
           blasDataType* y,
           const int* incY);

// SYR
void SYR_(const char* Uplo,
          const int* N,
          const blasDataType* alpha,
          const blasDataType* X,
          const int* incX,
          blasDataType* A,
          const int* lda);

// GER
void GER_(const int* M,
          const int* N,
          const blasDataType* alpha,
          const blasDataType* X,
          const int* incX,
          const blasDataType* Y,
          const int* incY,
          blasDataType* A,
          const int* lda);

// GEMV
void GEMV_(const char* trans,
           const int* m,
           const int* n,
           const blasDataType* alpha,
           const blasDataType* a,
           const int* lda,
           const blasDataType* x,
           const int* incX,
           const blasDataType* beta,
           blasDataType* y,
           const int* incY);

// GEMM
void GEMM_(const char* TransA,
           const char* TransB,
           const int* M,
           const int* N,
           const int* K,
           const blasDataType* alpha,
           const blasDataType* A,
           const int* lda,
           const blasDataType* B,
           const int* ldb,
           const blasDataType* beta,
           blasDataType* C,
           const int* ldc);
// TRSM
void TRSM_(const char* Side,
           const char* Uplo,
           const char* TransA,
           const char* Diag,
           const int* M,
           const int* N,
           const blasDataType* alpha,
           const blasDataType* A,
           const int* lda,
           blasDataType* B,
           const int* ldb);