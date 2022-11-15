/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTRS_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-27
 *******************************************************************************/

#pragma once
#include <complex.h>
#include <kblas.h>
#include <kml_service.h>
#include <math.h>
#include "type.h"

/*
void SYTRS_ROOK(const char* uplo,
                const int* n,
                const int* nrhs,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* B,
                const int* ldb,
                int* info);
*/

void ssytrs_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    float* A,
                    const int* lda,
                    int* ipiv,
                    float* B,
                    const int* ldb,
                    int* info);

void dsytrs_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    double* A,
                    const int* lda,
                    int* ipiv,
                    double* B,
                    const int* ldb,
                    int* info);

void csytrs_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    kml_complex_float* A,
                    const int* lda,
                    int* ipiv,
                    kml_complex_float* B,
                    const int* ldb,
                    int* info);

void zsytrs_rook_(const char* uplo,
                    const int* n,
                    const int* nrhs,
                    kml_complex_double* A,
                    const int* lda,
                    int* ipiv,
                    kml_complex_double* B,
                    const int* ldb,
                    int* info);
                    