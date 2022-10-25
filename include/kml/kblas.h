/*
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights reserved.
 * Description: cblas interface for the user
 * Author: KML
 * Create: 2020-09-22
 * Notes: NA
 */

#ifndef CBLAS_H
#define CBLAS_H

#include <stddef.h>

#ifndef BLAS_TYPES_DEF
#define BLAS_TYPES_DEF
#ifndef ASSEMBLER
typedef long BLASLONG;
typedef unsigned long BLASULONG;
#ifdef USE64BITINT
typedef BLASLONG BLASINT;
#define blasabs(x) labs(x)
#else
typedef int BLASINT;
#define blasabs(x) abs(x)
#endif // USE64BITINT

/* C99 supports complex floating numbers natively, which GCC also offers as an
   extension since version 3.0.  If neither are available, use a compatible
   structure as fallback (see Clause 6.2.5.13 of the C99 standard). */
#if ((defined(__STDC_IEC_559_COMPLEX__) || __STDC_VERSION__ >= 199901L || \
    (__GNUC__ >= 3 && !defined(__cplusplus))) && \
    !(defined(FORCE_BLAS_COMPLEX_STRUCT)))
#ifndef __cplusplus
#include <complex.h>
#endif
typedef float _Complex blas_complex_float;
typedef double _Complex blas_complex_double;
#define blas_make_complex_float(real, imag) ((real) + ((imag)*_Complex_I))
#define blas_make_complex_double(real, imag) ((real) + ((imag)*_Complex_I))
#else
#error("Error:required:__STDC_IEC_559_COMPLEX__ || __STDC_VERSION__ >= 199901L || __GNUC__ >= 3")
#endif
#endif // ASSEMBLER
#endif // BLAS_TYPES_DEF

#ifdef __cplusplus
extern "C" {
/* Assume C declarations for C++ */
#endif /* __cplusplus */

/* Set the number of threads on runtime. */
void BlasSetNumThreads(int numThreads);

/* Get the number of threads on runtime. */
int BlasGetNumThreads(void);

/* Get the number of physical processors (cores). */
int BlasGetNumProcs(void);

/* KML_BLAS is compiled for sequential use  */
#define BLAS_SEQUENTIAL 0
/* KML_BLAS is compiled using normal threading model */
#define BLAS_THREAD 1
/* KML_BLAS is compiled using OpenMP threading model */
#define BLAS_OPENMP 2

/* Get KML_BLAS parallelization type */
int BlasGetParallel(void);

#define CBLAS_INDEX size_t

/* matrix saved in rows or cols */
typedef enum CBLAS_ORDER {
    CblasRowMajor = 101,
    CblasColMajor = 102
} CBLAS_ORDER;

/* matrix transpose or conjugate transpose */
typedef enum CBLAS_TRANSPOSE {
    CblasNoTrans = 111,
    CblasTrans = 112,
    CblasConjTrans = 113, // conjugate transpose
    CblasConjNoTrans = 114
} CBLAS_TRANSPOSE;

/* matrix upper or lower */
typedef enum CBLAS_UPLO {
    CblasUpper = 121,
    CblasLower = 122
} CBLAS_UPLO;

/* matrix diagonal is unit or no unit */
typedef enum CBLAS_DIAG {
    CblasNonUnit = 131,
    CblasUnit = 132
} CBLAS_DIAG;

/* matrix position is left or right */
typedef enum CBLAS_SIDE {
    CblasLeft = 141,
    CblasRight = 142
} CBLAS_SIDE;
typedef CBLAS_ORDER CBLAS_LAYOUT;

#define VERSION_STRUCT_LEN 100
typedef struct {
    char component[VERSION_STRUCT_LEN];
    char version[VERSION_STRUCT_LEN];
    char supportPackage[VERSION_STRUCT_LEN];
    char compileTime[VERSION_STRUCT_LEN];
}KBLASVersion;
int KBLASGetVersion(KBLASVersion *ver);

float cblas_sdsdot(const BLASINT n, const float alpha, const float *x, const BLASINT incX, const float *y,
    const BLASINT incY);
double cblas_dsdot(const BLASINT n, const float *x, const BLASINT incX, const float *y, const BLASINT incY);
float cblas_sdot(const BLASINT n, const float *x, const BLASINT incX, const float *y, const BLASINT incY);
double cblas_ddot(const BLASINT n, const double *x, const BLASINT incX, const double *y, const BLASINT incY);

blas_complex_float cblas_cdotu(const BLASINT n, const void *x, const BLASINT incX, const void *y,
    const BLASINT incY);
blas_complex_float cblas_cdotc(const BLASINT n, const void *x, const BLASINT incX, const void *y,
    const BLASINT incY);
blas_complex_double cblas_zdotu(const BLASINT n, const void *x, const BLASINT incX, const void *y,
    const BLASINT incY);
blas_complex_double cblas_zdotc(const BLASINT n, const void *x, const BLASINT incX, const void *y,
    const BLASINT incY);

void cblas_cdotu_sub(const BLASINT n, const void *x, const BLASINT incX, const void *y, const BLASINT incY, void *ret);
void cblas_cdotc_sub(const BLASINT n, const void *x, const BLASINT incX, const void *y, const BLASINT incY, void *ret);
void cblas_zdotu_sub(const BLASINT n, const void *x, const BLASINT incX, const void *y, const BLASINT incY, void *ret);
void cblas_zdotc_sub(const BLASINT n, const void *x, const BLASINT incX, const void *y, const BLASINT incY, void *ret);

float cblas_sasum(const BLASINT n, const float *x, const BLASINT incX);
double cblas_dasum(const BLASINT n, const double *x, const BLASINT incX);
float cblas_scasum(const BLASINT n, const void *x, const BLASINT incX);
double cblas_dzasum(const BLASINT n, const void *x, const BLASINT incX);

float cblas_snrm2(const BLASINT N, const float *X, const BLASINT incX);
double cblas_dnrm2(const BLASINT N, const double *X, const BLASINT incX);
float cblas_scnrm2(const BLASINT N, const void *X, const BLASINT incX);
double cblas_dznrm2(const BLASINT N, const void *X, const BLASINT incX);

CBLAS_INDEX cblas_isamax(const BLASINT n, const float *x, const BLASINT incX);
CBLAS_INDEX cblas_idamax(const BLASINT n, const double *x, const BLASINT incX);
CBLAS_INDEX cblas_icamax(const BLASINT n, const void *x, const BLASINT incX);
CBLAS_INDEX cblas_izamax(const BLASINT n, const void *x, const BLASINT incX);

CBLAS_INDEX cblas_isamin(const BLASINT n, const float *x, const BLASINT incX);
CBLAS_INDEX cblas_idamin(const BLASINT n, const double *x, const BLASINT incX);
CBLAS_INDEX cblas_icamin(const BLASINT n, const void *x, const BLASINT incX);
CBLAS_INDEX cblas_izamin(const BLASINT n, const void *x, const BLASINT incX);

void cblas_saxpy(const BLASINT n, const float alpha, const float *x, const BLASINT incX, float *y, const BLASINT incY);
void cblas_daxpy(const BLASINT n, const double alpha, const double *x, const BLASINT incX, double *y,
    const BLASINT incY);
void cblas_caxpy(const BLASINT n, const void *alpha, const void *x, const BLASINT incX, void *y, const BLASINT incY);
void cblas_zaxpy(const BLASINT n, const void *alpha, const void *x, const BLASINT incX, void *y, const BLASINT incY);
void cblas_caxpyc(const BLASINT n, const void *alpha, const void *x, const BLASINT incX, void *y, const BLASINT incY);
void cblas_zaxpyc(const BLASINT n, const void *alpha, const void *x, const BLASINT incX, void *y, const BLASINT incY);

void cblas_scopy(const BLASINT n, const float *x, const BLASINT incX, float *y, const BLASINT incY);
void cblas_dcopy(const BLASINT n, const double *x, const BLASINT incX, double *y, const BLASINT incY);
void cblas_ccopy(const BLASINT n, const void *x, const BLASINT incX, void *y, const BLASINT incY);
void cblas_zcopy(const BLASINT n, const void *x, const BLASINT incX, void *y, const BLASINT incY);

void cblas_sswap(const BLASINT n, float *x, const BLASINT incX, float *y, const BLASINT incY);
void cblas_dswap(const BLASINT n, double *x, const BLASINT incX, double *y, const BLASINT incY);
void cblas_cswap(const BLASINT n, void *x, const BLASINT incX, void *y, const BLASINT incY);
void cblas_zswap(const BLASINT n, void *x, const BLASINT incX, void *y, const BLASINT incY);

void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_drotg(double *a, double *b, double *c, double *s);

void cblas_srot(const BLASINT N, float *X, const BLASINT incX, float *Y, const BLASINT incY, const float c,
    const float s);
void cblas_drot(const BLASINT N, double *X, const BLASINT incX, double *Y, const BLASINT incY, const double c,
    const double s);
void cblas_csrot(const BLASINT N, float *X, const BLASINT incX, float *Y, const BLASINT incY, const float c,
    const float s);
void cblas_zdrot(const BLASINT N, double *X, const BLASINT incX, double *Y, const BLASINT incY, const double c,
    const double s);

void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);

void cblas_srotm(const BLASINT N, float *X, const BLASINT incX, float *Y, const BLASINT incY, const float *P);
void cblas_drotm(const BLASINT N, double *X, const BLASINT incX, double *Y, const BLASINT incY, const double *P);

void cblas_sscal(const BLASINT N, const float alpha, float *X, const BLASINT incX);
void cblas_dscal(const BLASINT N, const double alpha, double *X, const BLASINT incX);
void cblas_cscal(const BLASINT N, const void *alpha, void *X, const BLASINT incX);
void cblas_zscal(const BLASINT N, const void *alpha, void *X, const BLASINT incX);
void cblas_csscal(const BLASINT N, const float alpha, void *X, const BLASINT incX);
void cblas_zdscal(const BLASINT N, const double alpha, void *X, const BLASINT incX);

void cblas_sgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const BLASINT m, const BLASINT n,
    const float alpha, const float *a, const BLASINT lda, const float *x, const BLASINT incX, const float beta,
    float *y, const BLASINT incY);
void cblas_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const BLASINT m, const BLASINT n,
    const double alpha, const double *a, const BLASINT lda, const double *x, const BLASINT incX, const double beta,
    double *y, const BLASINT incY);
void cblas_cgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const BLASINT m, const BLASINT n,
    const void *alpha, const void *a, const BLASINT lda, const void *x, const BLASINT incX, const void *beta, void *y,
    const BLASINT incY);
void cblas_zgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const BLASINT m, const BLASINT n,
    const void *alpha, const void *a, const BLASINT lda, const void *x, const BLASINT incX, const void *beta, void *y,
    const BLASINT incY);

void cblas_sger(const enum CBLAS_ORDER order, const BLASINT M, const BLASINT N, const float alpha, const float *X,
    const BLASINT incX, const float *Y, const BLASINT incY, float *A, const BLASINT lda);
void cblas_dger(const enum CBLAS_ORDER order, const BLASINT M, const BLASINT N, const double alpha, const double *X,
    const BLASINT incX, const double *Y, const BLASINT incY, double *A, const BLASINT lda);
void cblas_cgeru(const enum CBLAS_ORDER order, const BLASINT M, const BLASINT N, const void *alpha, const void *X,
    const BLASINT incX, const void *Y, const BLASINT incY, void *A, const BLASINT lda);
void cblas_cgerc(const enum CBLAS_ORDER order, const BLASINT M, const BLASINT N, const void *alpha, const void *X,
    const BLASINT incX, const void *Y, const BLASINT incY, void *A, const BLASINT lda);
void cblas_zgeru(const enum CBLAS_ORDER order, const BLASINT M, const BLASINT N, const void *alpha, const void *X,
    const BLASINT incX, const void *Y, const BLASINT incY, void *A, const BLASINT lda);
void cblas_zgerc(const enum CBLAS_ORDER order, const BLASINT M, const BLASINT N, const void *alpha, const void *X,
    const BLASINT incX, const void *Y, const BLASINT incY, void *A, const BLASINT lda);

void cblas_strsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const float *A, const BLASINT lda, float *X, const BLASINT incX);
void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const double *A, const BLASINT lda, double *X, const BLASINT incX);
void cblas_ctrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *A, const BLASINT lda, void *X, const BLASINT incX);
void cblas_ztrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *A, const BLASINT lda, void *X, const BLASINT incX);

void cblas_strmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const float *A, const BLASINT lda, float *X, const BLASINT incX);
void cblas_dtrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const double *A, const BLASINT lda, double *X, const BLASINT incX);
void cblas_ctrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *A, const BLASINT lda, void *X, const BLASINT incX);
void cblas_ztrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *A, const BLASINT lda, void *X, const BLASINT incX);

void cblas_ssyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const float *X, const BLASINT incX, float *A, const BLASINT lda);
void cblas_dsyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const double *X, const BLASINT incX, double *A, const BLASINT lda);
void cblas_cher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const void *X, const BLASINT incX, void *A, const BLASINT lda);
void cblas_zher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const void *X, const BLASINT incX, void *A, const BLASINT lda);

void cblas_ssyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const float *X, const BLASINT incX, const float *Y, const BLASINT incY, float *A, const BLASINT lda);
void cblas_dsyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const double *X, const BLASINT incX, const double *Y, const BLASINT incY, double *A, const BLASINT lda);
void cblas_cher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *X, const BLASINT incX, const void *Y, const BLASINT incY, void *A, const BLASINT lda);
void cblas_zher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *X, const BLASINT incX, const void *Y, const BLASINT incY, void *A, const BLASINT lda);

void cblas_sgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const BLASINT M, const BLASINT N,
    const BLASINT KL, const BLASINT KU, const float alpha, const float *A, const BLASINT lda, const float *X,
    const BLASINT incX, const float beta, float *Y, const BLASINT incY);
void cblas_dgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const BLASINT M, const BLASINT N,
    const BLASINT KL, const BLASINT KU, const double alpha, const double *A, const BLASINT lda, const double *X,
    const BLASINT incX, const double beta, double *Y, const BLASINT incY);
void cblas_cgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const BLASINT M, const BLASINT N,
    const BLASINT KL, const BLASINT KU, const void *alpha, const void *A, const BLASINT lda, const void *X,
    const BLASINT incX, const void *beta, void *Y, const BLASINT incY);
void cblas_zgbmv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA, const BLASINT M, const BLASINT N,
    const BLASINT KL, const BLASINT KU, const void *alpha, const void *A, const BLASINT lda, const void *X,
    const BLASINT incX, const void *beta, void *Y, const BLASINT incY);

void cblas_ssbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const BLASINT K,
    const float alpha, const float *A, const BLASINT lda, const float *X, const BLASINT incX, const float beta,
    float *Y, const BLASINT incY);
void cblas_dsbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const BLASINT K,
    const double alpha, const double *A, const BLASINT lda, const double *X, const BLASINT incX, const double beta,
    double *Y, const BLASINT incY);


void cblas_stbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const float *A, const BLASINT lda, float *X,
    const BLASINT incX);
void cblas_dtbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const double *A, const BLASINT lda, double *X,
    const BLASINT incX);
void cblas_ctbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const void *A, const BLASINT lda, void *X,
    const BLASINT incX);
void cblas_ztbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const void *A, const BLASINT lda, void *X,
    const BLASINT incX);

void cblas_stbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const float *A, const BLASINT lda, float *X,
    const BLASINT incX);
void cblas_dtbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const double *A, const BLASINT lda, double *X,
    const BLASINT incX);
void cblas_ctbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const void *A, const BLASINT lda, void *X,
    const BLASINT incX);
void cblas_ztbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const BLASINT K, const void *A, const BLASINT lda, void *X,
    const BLASINT incX);

void cblas_stpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const float *Ap, float *X, const BLASINT incX);
void cblas_dtpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const double *Ap, double *X, const BLASINT incX);
void cblas_ctpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *Ap, void *X, const BLASINT incX);
void cblas_ztpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *Ap, void *X, const BLASINT incX);

void cblas_stpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const float *Ap, float *X, const BLASINT incX);
void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const double *Ap, double *X, const BLASINT incX);
void cblas_ctpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *Ap, void *X, const BLASINT incX);
void cblas_ztpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag, const BLASINT N, const void *Ap, void *X, const BLASINT incX);

void cblas_ssymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const float *A, const BLASINT lda, const float *X, const BLASINT incX, const float beta, float *Y,
    const BLASINT incY);
void cblas_dsymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const double *A, const BLASINT lda, const double *X, const BLASINT incX, const double beta, double *Y,
    const BLASINT incY);
void cblas_chemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *A, const BLASINT lda, const void *X, const BLASINT incX, const void *beta, void *Y, const BLASINT incY);
void cblas_zhemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *A, const BLASINT lda, const void *X, const BLASINT incX, const void *beta, void *Y, const BLASINT incY);


void cblas_sspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const float *Ap, const float *X, const BLASINT incX, const float beta, float *Y, const BLASINT incY);
void cblas_dspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const double *Ap, const double *X, const BLASINT incX, const double beta, double *Y, const BLASINT incY);

void cblas_sspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const float *X, const BLASINT incX, float *Ap);
void cblas_dspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const double *X, const BLASINT incX, double *Ap);

void cblas_chpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const void *X, const BLASINT incX, void *A);
void cblas_zhpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const void *X, const BLASINT incX, void *A);

void cblas_sspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const float alpha,
    const float *X, const BLASINT incX, const float *Y, const BLASINT incY, float *A);
void cblas_dspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const double alpha,
    const double *X, const BLASINT incX, const double *Y, const BLASINT incY, double *A);
void cblas_chpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *X, const BLASINT incX, const void *Y, const BLASINT incY, void *Ap);
void cblas_zhpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *X, const BLASINT incX, const void *Y, const BLASINT incY, void *Ap);

void cblas_chbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const BLASINT K,
    const void *alpha, const void *A, const BLASINT lda, const void *X, const BLASINT incX, const void *beta, void *Y,
    const BLASINT incY);
void cblas_zhbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const BLASINT K,
    const void *alpha, const void *A, const BLASINT lda, const void *X, const BLASINT incX, const void *beta, void *Y,
    const BLASINT incY);

void cblas_chpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *Ap, const void *X, const BLASINT incX, const void *beta, void *Y, const BLASINT incY);
void cblas_zhpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const BLASINT N, const void *alpha,
    const void *Ap, const void *X, const BLASINT incX, const void *beta, void *Y, const BLASINT incY);

void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
    const BLASINT M, const BLASINT N, const BLASINT K, const float alpha, const float *A, const BLASINT lda,
    const float *B, const BLASINT ldb, const float beta, float *C, const BLASINT ldc);
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
    const BLASINT M, const BLASINT N, const BLASINT K, const double alpha, const double *A, const BLASINT lda,
    const double *B, const BLASINT ldb, const double beta, double *C, const BLASINT ldc);
void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
    const BLASINT M, const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda,
    const void *B, const BLASINT ldb, const void *beta, void *C, const BLASINT ldc);
void cblas_cgemm3m(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
    const BLASINT M, const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda,
    const void *B, const BLASINT ldb, const void *beta, void *C, const BLASINT ldc);
void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
    const BLASINT M, const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda,
    const void *B, const BLASINT ldb, const void *beta, void *C, const BLASINT ldc);
void cblas_zgemm3m(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB,
    const BLASINT M, const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda,
    const void *B, const BLASINT ldb, const void *beta, void *C, const BLASINT ldc);


void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const BLASINT M,
    const BLASINT N, const float alpha, const float *A, const BLASINT lda, const float *B, const BLASINT ldb,
    const float beta, float *C, const BLASINT ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const BLASINT M,
    const BLASINT N, const double alpha, const double *A, const BLASINT lda, const double *B, const BLASINT ldb,
    const double beta, double *C, const BLASINT ldc);
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const BLASINT M,
    const BLASINT N, const void *alpha, const void *A, const BLASINT lda, const void *B, const BLASINT ldb,
    const void *beta, void *C, const BLASINT ldc);
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const BLASINT M,
    const BLASINT N, const void *alpha, const void *A, const BLASINT lda, const void *B, const BLASINT ldb,
    const void *beta, void *C, const BLASINT ldc);

void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const float alpha, const float *A, const BLASINT lda, const float beta, float *C,
    const BLASINT ldc);
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const double alpha, const double *A, const BLASINT lda, const double beta,
    double *C, const BLASINT ldc);
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda, const void *beta, void *C,
    const BLASINT ldc);
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda, const void *beta, void *C,
    const BLASINT ldc);

void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const float alpha, const float *A, const BLASINT lda, const float *B,
    const BLASINT ldb, const float beta, float *C, const BLASINT ldc);
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const double alpha, const double *A, const BLASINT lda, const double *B,
    const BLASINT ldb, const double beta, double *C, const BLASINT ldc);
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda, const void *B,
    const BLASINT ldb, const void *beta, void *C, const BLASINT ldc);
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda, const void *B,
    const BLASINT ldb, const void *beta, void *C, const BLASINT ldc);

void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const float alpha,
    const float *A, const BLASINT lda, float *B, const BLASINT ldb);
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const double alpha,
    const double *A, const BLASINT lda, double *B, const BLASINT ldb);
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const void *alpha,
    const void *A, const BLASINT lda, void *B, const BLASINT ldb);
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const void *alpha,
    const void *A, const BLASINT lda, void *B, const BLASINT ldb);

void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const float alpha,
    const float *A, const BLASINT lda, float *B, const BLASINT ldb);
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const double alpha,
    const double *A, const BLASINT lda, double *B, const BLASINT ldb);
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const void *alpha,
    const void *A, const BLASINT lda, void *B, const BLASINT ldb);
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag, const BLASINT M, const BLASINT N, const void *alpha,
    const void *A, const BLASINT lda, void *B, const BLASINT ldb);

void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const BLASINT M,
    const BLASINT N, const void *alpha, const void *A, const BLASINT lda, const void *B, const BLASINT ldb,
    const void *beta, void *C, const BLASINT ldc);
void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side, const enum CBLAS_UPLO Uplo, const BLASINT M,
    const BLASINT N, const void *alpha, const void *A, const BLASINT lda, const void *B, const BLASINT ldb,
    const void *beta, void *C, const BLASINT ldc);

void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const float alpha, const void *A, const BLASINT lda, const float beta, void *C,
    const BLASINT ldc);
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const double alpha, const void *A, const BLASINT lda, const double beta, void *C,
    const BLASINT ldc);

void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda, const void *B,
    const BLASINT ldb, const float beta, void *C, const BLASINT ldc);
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE Trans,
    const BLASINT N, const BLASINT K, const void *alpha, const void *A, const BLASINT lda, const void *B,
    const BLASINT ldb, const double beta, void *C, const BLASINT ldc);

void cblas_xerbla(BLASINT p, const char *rout, const char *form, ...);

/* ** BLAS extensions ** */

void cblas_saxpby(const BLASINT n, const float alpha, const float *x, const BLASINT incX, const float beta, float *y,
    const BLASINT incY);

void cblas_daxpby(const BLASINT n, const double alpha, const double *x, const BLASINT incX, const double beta,
    double *y, const BLASINT incY);

void cblas_caxpby(const BLASINT n, const void *alpha, const void *x, const BLASINT incX, const void *beta, void *y,
    const BLASINT incY);

void cblas_zaxpby(const BLASINT n, const void *alpha, const void *x, const BLASINT incX, const void *beta, void *y,
    const BLASINT incY);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif
