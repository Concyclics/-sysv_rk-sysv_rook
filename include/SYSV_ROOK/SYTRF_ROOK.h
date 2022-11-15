/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTRF_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-29
 *******************************************************************************/

#pragma once

#include <complex.h>
#include <kblas.h>
#include <kml_service.h>
#include <math.h>
#include "LASYF_ROOK.h"
#include "SYTF2_ROOK.h"
#include "ilaenv.h"
#include "type.h"

/*
void SYTRF_ROOK(const char* uplo,
                const int* n,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* work,
                const int* lwork,
                int* info);
*/

void ssytrf_rook_(const char* uplo,
                    const int* n,
                    float* A,
                    const int* lda,
                    int* ipiv,
                    float* work,
                    const int* lwork,
                    int* info);

void dsytrf_rook_(const char* uplo,
                    const int* n,
                    double* A,
                    const int* lda,
                    int* ipiv,
                    double* work,
                    const int* lwork,
                    int* info);

void csytrf_rook_(const char* uplo,
                    const int* n,
                    kml_complex_float* A,
                    const int* lda,
                    int* ipiv,
                    kml_complex_float* work,
                    const int* lwork,
                    int* info);

void zsytrf_rook_(const char* uplo,
                    const int* n,
                    kml_complex_double* A,
                    const int* lda,
                    int* ipiv,
                    kml_complex_double* work,
                    const int* lwork,
                    int* info);
                    