/*
* Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights reserved.
* Description: KML function
* Author: KML
* Create: 2020
*/
#ifndef KLAPACK_EXTERNAL_H
#define KLAPACK_EXTERNAL_H

#include <complex.h>
#include <kml_export.h>

KML_EXPORT void cppsv_(
    const char *uplo,
    const int *n, const int *nrhs,
    float _Complex *ap,
    float _Complex *b, const int *ldb,
    int *info);

KML_EXPORT void dppsv_(
    const char *uplo,
    const int *n, const int *nrhs,
    double *ap,
    double *b, const int *ldb,
    int *info);

KML_EXPORT void sppsv_(
    const char *uplo,
    const int *n, const int *nrhs,
    float *ap,
    float *b, const int *ldb,
    int *info);

KML_EXPORT void zppsv_(
    const char *uplo,
    const int *n, const int *nrhs,
    double _Complex *ap,
    double _Complex *b, const int *ldb,
    int *info);

KML_EXPORT void cgetrf_(
    const int *m, const int *n,
    float _Complex *a, const int *lda, int *ipiv,
    int *info);

KML_EXPORT void dgetrf_(
    const int *m, const int *n,
    double *a, const int *lda, int *ipiv,
    int *info);

KML_EXPORT void sgetrf_(
    const int *m, const int *n,
    float *a, const int *lda, int *ipiv,
    int *info);

KML_EXPORT void zgetrf_(
    const int *m, const int *n,
    double _Complex *a, const int *lda, int *ipiv,
    int *info);

KML_EXPORT void cgetri_(
    const int *n,
    float _Complex *a, const int *lda, const int *ipiv,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void dgetri_(
    const int *n,
    double *a, const int *lda, const int *ipiv,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void sgetri_(
    const int *n,
    float *a, const int *lda, const int *ipiv,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void zgetri_(
    const int *n,
    double _Complex *a, const int *lda, const int *ipiv,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void cgeqrf_(
    const int *m, const int *n,
    float _Complex *a, const int *lda,
    float _Complex *tau,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void dgeqrf_(
    const int *m, const int *n,
    double *a, const int *lda,
    double *tau,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void sgeqrf_(
    const int *m, const int *n,
    float *a, const int *lda,
    float *tau,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void zgeqrf_(
    const int *m, const int *n,
    double _Complex *a, const int *lda,
    double _Complex *tau,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void sorglq_(
    const int *m, const int *n, const int *k,
    float *a, const int *lda,
    const float *tau,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dorglq_(
    const int *m, const int *n, const int *k,
    double *a, const int *lda,
    const double *tau,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void cunglq_(
    const int *m, const int *n, const int *k,
    float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zunglq_(
    const int *m, const int *n, const int *k,
    double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void sorgql_(
    const int *m, const int *n, const int *k,
    float *a, const int *lda,
    const float *tau,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dorgql_(
    const int *m, const int *n, const int *k,
    double *a, const int *lda,
    const double *tau,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void cungql_(
    const int *m, const int *n, const int *k,
    float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zungql_(
    const int *m, const int *n, const int *k,
    double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void sorgqr_(
    const int *m, const int *n, const int *k,
    float *a, const int *lda,
    const float *tau,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dorgqr_(
    const int *m, const int *n, const int *k,
    double *a, const int *lda,
    const double *tau,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void cungqr_(
    const int *m, const int *n, const int *k,
    float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zungqr_(
    const int *m, const int *n, const int *k,
    double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void sorgrq_(
    const int *m, const int *n, const int *k,
    float *a, const int *lda,
    const float *tau,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dorgrq_(
    const int *m, const int *n, const int *k,
    double *a, const int *lda,
    const double *tau,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void cungrq_(
    const int *m, const int *n, const int *k,
    float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zungrq_(
    const int *m, const int *n, const int *k,
    double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void dormqr_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double *a, const int *lda,
    const double *tau,
    double *c, const int *ldc,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void sormqr_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float *a, const int *lda,
    const float *tau,
    float *c, const int *ldc,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dormlq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double *a, const int *lda,
    const double *tau,
    double *c, const int *ldc,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void sormlq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float *a, const int *lda,
    const float *tau,
    float *c, const int *ldc,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dormql_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double *a, const int *lda,
    const double *tau,
    double *c, const int *ldc,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void sormql_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float *a, const int *lda,
    const float *tau,
    float *c, const int *ldc,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void dormrq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double *a, const int *lda,
    const double *tau,
    double *c, const int *ldc,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void sormrq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float *a, const int *lda,
    const float *tau,
    float *c, const int *ldc,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void cunmqr_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *c, const int *ldc,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zunmqr_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *c, const int *ldc,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void cunmlq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *c, const int *ldc,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zunmlq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *c, const int *ldc,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void cunmql_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *c, const int *ldc,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zunmql_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *c, const int *ldc,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void cunmrq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const float _Complex *a, const int *lda,
    const float _Complex *tau,
    float _Complex *c, const int *ldc,
    float _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void zunmrq_(
    const char *side, const char *trans,
    const int *m, const int *n, const int *k,
    const double _Complex *a, const int *lda,
    const double _Complex *tau,
    double _Complex *c, const int *ldc,
    double _Complex *work, const int *lwork,
    int *info);

KML_EXPORT void cpotrf_(
    const char *uplo,
    const int *n,
    float _Complex *a, const int *lda,
    int *info);

KML_EXPORT void dpotrf_(
    const char *uplo,
    const int *n,
    double *a, const int *lda,
    int *info);

KML_EXPORT void spotrf_(
    const char *uplo,
    const int *n,
    float *a, const int *lda,
    int *info);

KML_EXPORT void zpotrf_(
    const char *uplo,
    const int *n,
    double _Complex *a, const int *lda,
    int *info);

KML_EXPORT void cgesv_(
    const int *n, const int *nrhs,
    float _Complex *a, const int *lda, int *ipiv,
    float _Complex *b, const int *ldb,
    int *info);

KML_EXPORT void dgesv_(
    const int *n, const int *nrhs,
    double *a, const int *lda, int *ipiv,
    double *b, const int *ldb,
    int *info);

KML_EXPORT void sgesv_(
    const int *n, const int *nrhs,
    float *a, const int *lda, int *ipiv,
    float *b, const int *ldb,
    int *info);

KML_EXPORT void zgesv_(
    const int *n, const int *nrhs,
    double _Complex *a, const int *lda, int *ipiv,
    double _Complex *b, const int *ldb,
    int *info);

KML_EXPORT void ssyevd_(
    const char *jobz, const char *uplo,
    const int *n, float *a, const int *lda,
    float *w, float *work, const int *lwork,
    int *iwork, const int *liwork, int *info);

KML_EXPORT void dsyevd_(
    const char *jobz, const char *uplo,
    const int *n, double *a, const int *lda,
    double *w, double *work, const int *lwork,
    int *iwork, const int *liwork, int *info);

KML_EXPORT void cheevd_(
    const char *jobz, const char *uplo,
    const int *n, float _Complex *a, const int *lda,
    float *w, float _Complex *work, const int *lwork,
    float *rwork, const int *lrwork, int *iwork, 
    const int *liwork, int *info);

KML_EXPORT void zheevd_(
    const char *jobz, const char *uplo,
    const int *n, double _Complex *a, const int *lda,
    double *w, double _Complex *work, const int *lwork,
    double *rwork, const int *lrwork, int *iwork, 
    const int *liwork, int *info);

KML_EXPORT void cpotri_(
    const char *uplo,
    const int *n,
    float _Complex *a, const int *lda,
    int *info);

KML_EXPORT void dpotri_(
    const char *uplo,
    const int *n,
    double *a, const int *lda,
    int *info);

KML_EXPORT void spotri_(
    const char *uplo,
    const int *n,
    float *a, const int *lda,
    int *info);

KML_EXPORT void zpotri_(
    const char *uplo,
    const int *n,
    double _Complex *a, const int *lda,
    int *info);

KML_EXPORT void chetrd_(
    const char *uplo, const int *n, float _Complex *a, const int *lda,
    float *d, float *e, float _Complex *tau,
    float _Complex *work, const int *lwork, int *info);

KML_EXPORT void dsytrd_(
    const char *uplo, const int *n, double *a, const int *lda,
    double *d, double *e, double *tau, double *work,
    const int *lwork, int *info);

KML_EXPORT void ssytrd_(
    const char *uplo, const int *n, float *a, const int *lda,
    float *d, float *e, float *tau, float *work,
    const int *lwork, int *info);

KML_EXPORT void zhetrd_(
    const char *uplo, const int *n, double _Complex *a, const int *lda,
    double *d, double *e, double _Complex *tau,
    double _Complex *work, const int *lwork, int *info);

KML_EXPORT void dsyev_(
    const char *jobz, const char *uplo,
    const int *n,
    double *a, const int *lda,
    double *w,
    double *work, const int *lwork,
    int *info);

KML_EXPORT void ssyev_(
    const char *jobz, const char *uplo,
    const int *n,
    float *a, const int *lda,
    float *w,
    float *work, const int *lwork,
    int *info);

KML_EXPORT void cheev_(
    const char *jobz, const char *uplo,
    const int *n,
    float _Complex *a, const int *lda,
    float *w,
    float _Complex *work, const int *lwork,
    float *rwork,
    int *info);

KML_EXPORT void zheev_(
    const char *jobz, const char *uplo,
    const int *n,
    double _Complex *a, const int *lda,
    double *w,
    double _Complex *work, const int *lwork,
    double *rwork,
    int *info);

typedef struct KMLVersion KLAPACKVersion;
KML_EXPORT int KLAPACKGetVersion(KLAPACKVersion *ver);

#endif /* KLAPACK_EXTERNAL_H */
