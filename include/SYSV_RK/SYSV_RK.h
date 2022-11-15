/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYSV_RK.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

#include <kml_service.h>
#include "SYTRF_RK.h"
#include "SYTRS_3.h"
#include "type.h"

/*
void SYSV_RK(const char* uplo,
             const int* n,
             const int* nrhs,
             dataType* a,
             const int* lda,
             dataType* e,
             int* ipiv,
             dataType* b,
             const int* ldb,
             dataType* work,
             const int* lwork,
             int* info);
*/

void ssysv_rk_(const char* uplo,
                const int* n,
                const int* nrhs,
                float* a,
                const int* lda,
                float* e,
                int* ipiv,
                float* b,
                const int* ldb,
                float* work,
                const int* lwork,
                int* info);

void dsysv_rk_(const char* uplo,
                const int* n,
                const int* nrhs,
                double* a,
                const int* lda,
                double* e,
                int* ipiv,
                double* b,
                const int* ldb,
                double* work,
                const int* lwork,
                int* info);

void csysv_rk_(const char* uplo,
                const int* n,
                const int* nrhs,
                kml_complex_float* a,
                const int* lda,
                kml_complex_float* e,
                int* ipiv,
                kml_complex_float* b,
                const int* ldb,
                kml_complex_float* work,
                const int* lwork,
                int* info);

void zsysv_rk_(const char* uplo,
                const int* n,
                const int* nrhs,
                kml_complex_double* a,
                const int* lda,
                kml_complex_double* e,
                int* ipiv,
                kml_complex_double* b,
                const int* ldb,
                kml_complex_double* work,
                const int* lwork,
                int* info);