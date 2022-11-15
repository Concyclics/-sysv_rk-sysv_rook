/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTF2_ROOK.
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
void SYTF2_ROOK(const char* uplo,
                const int* n,
                dataType* A,
                const int* lda,
                int* ipiv,
                int* info);
*/

void ssytf2_rook_(const char* uplo,
                    const int* n,
                    float* A,
                    const int* lda,
                    int* ipiv,
                    int* info);

void dsytf2_rook_(const char* uplo,
                    const int* n,
                    double* A,
                    const int* lda,
                    int* ipiv,
                    int* info);

void csytf2_rook_(const char* uplo,
                    const int* n,
                    kml_complex_float* A,
                    const int* lda,
                    int* ipiv,
                    int* info);

void zsytf2_rook_(const char* uplo,
                    const int* n,
                    kml_complex_double* A,
                    const int* lda,
                    int* ipiv,
                    int* info);