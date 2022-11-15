/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYSV_ROOK.
 * Author: linhouzhong
 * Create: 2022-06-30
 *******************************************************************************/

#pragma once
#include <kml_service.h>
#include "SYTRF_ROOK.h"
#include "SYTRS_ROOK.h"
#include "type.h"

/*
void SYSV_ROOK(const char* uplo,
               const int* n,
               const int* nrhs,
               dataType* a,
               const int* lda,
               int* ipiv,
               dataType* b,
               const int* ldb,
               dataType* work,
               const int* lwork,
               int* info);
*/

void ssysv_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    float* a,
                    const int* lda,
                    int* ipiv,
                    float* b,
                    const int* ldb,
                    float* work,
                    const int* lwork,
                    int* info);

void dsysv_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    double* a,
                    const int* lda,
                    int* ipiv,
                    double* b,
                    const int* ldb,
                    double* work,
                    const int* lwork,
                    int* info);

void csysv_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    kml_complex_float* a,
                    const int* lda,
                    int* ipiv,
                    kml_complex_float* b,
                    const int* ldb,
                    kml_complex_float* work,
                    const int* lwork,
                    int* info);

void zsysv_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    kml_complex_double* a,
                    const int* lda,
                    int* ipiv,
                    kml_complex_double* b,
                    const int* ldb,
                    kml_complex_double* work,
                    const int* lwork,
                    int* info);
                    