/*
 * Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights reserved.
 * Description: DesignWare kunpeng libmath interface definition.
 * Author:
 * Create: 2020-08-28
 */

#ifndef KML_VML_H
#define KML_VML_H

#ifdef __cplusplus
extern "C" {
#endif

#define KVML_VERSION_STRUCT_LEN 100
typedef struct {
    char component[KVML_VERSION_STRUCT_LEN];
    char version[KVML_VERSION_STRUCT_LEN];
    char supportPackage[KVML_VERSION_STRUCT_LEN];
    char compiler[KVML_VERSION_STRUCT_LEN];
    char compileTime[KVML_VERSION_STRUCT_LEN];
}KVMLVersion;
int KVMLGetVersion(KVMLVersion* ver);

/**
 *  KML_VML FUNCTION ACCURACY CONTROL
 *  KML_HA - when KML_HA is set, high accuracy VML functions are called
 *  KML_LA - when KML_LA is set, low accuracy VML functions are called
 *  KML_EP - when KML_EP is set, enhanced performance VML functions are called
 * */
#define KML_LA 0x1
#define KML_HA 0x2
#define KML_EP 0x3

/**
 * @Brief Adds the elements of two vectors.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src1, src2 	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vsadd(const int len, const float *src1, const float *src2, float *dst);
void vdadd(const int len, const double *src1, const double *src2, double *dst);

/**
 * @Brief Subs the elements of two vectors.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src1, src2 	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vssub(const int len, const float *src1, const float *src2, float *dst);
void vdsub(const int len, const double *src1, const double *src2, double *dst);

/**
 * @Brief performs element by element squaring of the vector.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vssqr(const int len, const float *src, float *dst);
void vdsqr(const int len, const double *src, double *dst);

/**
 * @BriefComputes a square root of vector elements.
 * @param[in]       len         Number of elements in the vector
 * @param[in]       src         Pointers to the source vectors.
 * @param[out]      dst         Pointer to the destination vector.
 **/
void vssqrt(const int len, const float *src, float *dst);
void vdsqrt(const int len, const double *src, double *dst);

/**
 * @Brief Performs element by element multiplication of vector.
 * @param[in]       len         Number of elements in the vector
 * @param[in]       src1, src2  Pointers to the source vectors.
 * @param[out]      p _dst      Pointer to the destination vector.
 **/
void vsmul(const int len, const float *src1, const float *src2, float *dst);
void vdmul(const int len, const double *src1, const double *src2, double *dst);

/**
 * @Brief performs element by element division of vector.
 * @param[in]       len         Number of elements in the vector
 * @param[in]       src1, src2  Pointers to the source vectors.
 * @param[out]      p _dst      Pointer to the destination vector.
 **/
void vsdiv(const int len, const float *src1, const float *src2, float *dst);
void vddiv(const int len, const double *src1, const double *src2, double *dst);

/**
 * @Brief computes an exponential of vector elements.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vsexp(const int len, const float *src, float *dst);
void vdexp(const int len, const double *src, double *dst);

/**
 * @Brief computes ln of vector elements.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vsln(const int len, const float *src, float *dst);
void vdln(const int len, const double *src, double *dst);

/**
 * @Brief computes log10 of vector elements.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vslog10(const int len, const float *src, float *dst);
void vdlog10(const int len, const double *src, double *dst);

/**
 * @Computes tangent of vector elements.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vstan(const int len, const float *src, float *dst);
void vdtan(const int len, const double *src, double *dst);

/**
 * @Computes inverse tangent of vector elements.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vsatan(const int len, const float *src, float *dst);
void vdatan(const int len, const double *src, double *dst);

/**
 * @Computes inverse tangent of vector elements.
 * @param[in]		len			Number of elements in the vector
 * @param[in]  		src       	Pointers to the source vectors.
 * @param[out] 		dst        	Pointer to the destination vector.
 * */
void vsatan2(const int len, const float *src1, const float *src2, float *dst);
void vdatan2(const int len, const double *src1, const double *src2, double *dst);

/**
 * @Brief computes pow of vector elements.
 * @param[in]       len         Number of elements in the vector
 * @param[in]       src         Pointers to the source vectors.
 * @param[out]      dst         Pointer to the destination vector.
 * */
void vspow(const int len, const float *src1, const float *src2, float *dst);
void vdpow(const int len, const double *src1, const double *src2, double *dst);

/**
 * @Brief computes sine of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vssin(const int len, const float *src, float *dst);
void vdsin(const int len, const double *src, double *dst);

/**
 * @Brief computes cos of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vscos(const int len, const float *src, float *dst);
void vdcos(const int len, const double *src, double *dst);

/**
 * @Brief computes inv of vector elements.
 * @param[in]   	len         Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vsinv(const int len, const float *src, float *dst);
void vdinv(const int len, const double *src, double *dst);

/**
 * @Brief computes cos sin of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src     	Pointer to the source vector.
 * @param[out]  	sindst   	Pointer to the destination vector.
 * @param[out]  	cosdst   	Pointer to the destination vector.
 * */
void vssincos(const int len, const float *src, float *sindst, float *cosdst);
void vdsincos(const int len, const double *src, double *sindst, double *cosdst);

/**
 * @Brief computes sinh of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vssinh(const int len, const float *src, float *dst);
void vdsinh(const int len, const double *src, double *dst);

/**
 * @Brief computes cos of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vscosh(const int len, const float *src, float *dst);
void vdcosh(const int len, const double *src, double *dst);

/**
 * @Brief computes sinh of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vsasinh(const int len, const float *src, float *dst);
void vdasinh(const int len, const double *src, double *dst);

/**
 * @Brief computes cos of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vsacosh(const int len, const float *src, float *dst);
void vdacosh(const int len, const double *src, double *dst);

/**
 * @Brief computes tanh of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vstanh(const int len, const float *src, float *dst);
void vdtanh(const int len, const double *src, double *dst);

/**
 * @Brief computes atanh of vector elements.
 * @param[in]   	len    		Number of elements in the vector
 * @param[in]   	src   	    Pointers to the source vectors.
 * @param[out]  	dst   	    Pointer to the destination vector.
 * */
void vsatanh(const int len, const float *src, float *dst);
void vdatanh(const int len, const double *src, double *dst);

#ifdef __cplusplus
}
#endif

#endif