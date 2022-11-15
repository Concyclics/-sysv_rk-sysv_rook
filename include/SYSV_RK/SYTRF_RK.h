/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTRF_RK.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/
#pragma once

#include <complex.h>
#include <kblas.h>
#include <kml_service.h>
#include <math.h>
#include "LASYF_RK.h"
#include "SYTF2_RK.h"
#include "ilaenv.h"
#include "type.h"

/*
void SYTRF_RK(const char* uplo,
              const int* n,
              dataType* A,
              const int* lda,
              dataType* E,
              int* ipiv,
              dataType* work,
              const int* lwork,
              int* info);
*/

void ssytrf_rk_(const char* uplo,
                const int* n,
                float* A,
                const int* lda,
                float* E,
                int* ipiv,
                float* work,
                const int* lwork,
                int* info);

void dsytrf_rk_(const char* uplo,
                const int* n,
                double* A,
                const int* lda,
                double* E,
                int* ipiv,
                double* work,
                const int* lwork,
                int* info);

void csytrf_rk_(const char* uplo,
                const int* n,
                kml_complex_float* A,
                const int* lda,
                kml_complex_float* E,
                int* ipiv,
                kml_complex_float* work,
                const int* lwork,
                int* info);

void zsytrf_rk_(const char* uplo,
                const int* n,
                kml_complex_double* A,
                const int* lda,
                kml_complex_double* E,
                int* ipiv,
                kml_complex_double* work,
                const int* lwork,
                int* info);