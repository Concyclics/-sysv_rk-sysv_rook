/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTF2_RK.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

#include <complex.h>
#include <kblas.h>
#include <kml_service.h>
#include <math.h>
#include "type.h"

/*
void SYTF2_RK(const char* uplo,
              const int* n,
              dataType* A,
              const int* lda,
              dataType* E,
              int* ipiv,
              int* info);
*/

void ssytf2_rk_(const char* uplo,
                const int* n,
                float* A,
                const int* lda,
                float* E,
                int* ipiv,
                int* info);

void dsytf2_rk_(const char* uplo,
                const int* n,
                double* A,
                const int* lda,
                double* E,
                int* ipiv,
                int* info);

void csytf2_rk_(const char* uplo,
                const int* n,
                kml_complex_float* A,
                const int* lda,
                kml_complex_float* E,
                int* ipiv,
                int* info);

void zsytf2_rk_(const char* uplo,
                const int* n,
                kml_complex_double* A,
                const int* lda,
                kml_complex_double* E,
                int* ipiv,
                int* info);
                