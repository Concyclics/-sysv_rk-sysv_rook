/*
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights reserved.
 * Description: DesignWare kunpeng libmath interface definition.
 * Author:
 * Create: 2020-08-28
 */

#ifndef KML_LIBM_H
#define KML_LIBM_H

#ifdef __cplusplus
extern "C" {
#endif

#define complex _Complex
#define _Complex_I  (__extension__  1.0iF)
#undef I
#define I _Complex_I

#define KM_VERSION_STRUCT_LEN 100
typedef struct {
    char component[KM_VERSION_STRUCT_LEN];
    char version[KM_VERSION_STRUCT_LEN];
    char supportPackage[KM_VERSION_STRUCT_LEN];
    char compiler[KM_VERSION_STRUCT_LEN];
    char compileTime[KM_VERSION_STRUCT_LEN];
}KMVersion;

int KMGetVersion(KMVersion* ver);
int KML9GetVersion(KMVersion* ver);
float crealf(float complex);
float cimagf(float complex);
float complex catanf(float complex);
float complex csinf(float complex);
float complex ctanf(float complex);
float complex ccosf(float complex);

double creal(double complex);
double cimag(double complex);
double complex  catan (double complex);
double complex  ccos (double complex);
double complex  csin (double complex);
double complex  ctan (double complex);

float asindf(float x);
float cosdf(float x);
float tgammaf(float x);
float expf(float x);
float exp2f(float x);
float logf(float x);
float log2f(float x);
float log10f(float x);
float powf(float x, float y);
float sinf(float x);
float cosf(float x);
float tanf(float x);
float atanf(float x);
float atan2f(float x, float y);
float sqrtf(float x);
float cbrtf(float x);
float sinhf(float x);
float coshf(float x);
float tanhf(float x);
float asinhf(float x);
float acoshf(float x);
void sincosf(float x, float *sinp, float *cosp);
float acosf(float x);
float asinf(float x);
float atanhf(float x);
float fmodf(float x, float y);

double exp(double x);
double exp2(double x);
double log(double x);
double log2(double x);
double log10(double x);
double pow(double x, double y);
double sin(double x);
double cos(double x);
double tan(double x);
double atan(double x);
double atan2(double y, double x);
double sqrt(double x);
double cbrt(double x);
double sinh(double x);
double cosh(double x);
double tanh(double x);
double asinh(double x);
double acosh(double x);
void sincos(double x, double *sinv, double *cosv);
double acos(double x);
double asin(double x);
double atanh(double x);

float tgammaf_18(float x);
double cbrt_18(double x);
double log_18(double x);
double log10_18(double x);

#ifdef __cplusplus
}
#endif

#endif
