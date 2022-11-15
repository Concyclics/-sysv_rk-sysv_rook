/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTRS_3.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

#include <complex.h>
#include <kblas.h>
#include <kml_service.h>
#include <math.h>
#include "ilaenv.h"
#include "type.h"

/*
void SYTRS_3(const char* uplo,
             const int* n,
             const int* nrhs,
             dataType* A,
             const int* lda,
             dataType* E,
             int* ipiv,
             dataType* B,
             const int* ldb,
             int* info);
*/

void ssytrs_3_(const char* uplo,
                const int* n,
                const int* nrhs,
                float* A,
                const int* lda,
                float* E,
                int* ipiv,
                float* B,
                const int* ldb,
                int* info);

void dsytrs_3_(const char* uplo,
                const int* n,
                const int* nrhs,
                double* A,
                const int* lda,
                double* E,
                int* ipiv,
                double* B,
                const int* ldb,
                int* info);

void csytrs_3_(const char* uplo,
                const int* n,
                const int* nrhs,
                kml_complex_float* A,
                const int* lda,
                kml_complex_float* E,
                int* ipiv,
                kml_complex_float* B,
                const int* ldb,
                int* info);

void zsytrs_3_(const char* uplo,
                const int* n,
                const int* nrhs,
                kml_complex_double* A,
                const int* lda,
                kml_complex_double* E,
                int* ipiv,
                kml_complex_double* B,
                const int* ldb,
                int* info);
                