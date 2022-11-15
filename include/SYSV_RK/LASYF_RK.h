/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of LASYF_RK.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once
#include <complex.h>
#include <kblas.h>
#include <kml_service.h>
#include <math.h>
#include "type.h"

/*
void LASYF_RK(const char* uplo,
              const int* n,
              const int* nb,
              int* kb,
              dataType* a,
              const int* lda,
              dataType* e,
              int* ipiv,
              dataType* w,
              const int* ldw,
              int* info);
*/

void slasyf_rk_(const char* uplo,
                const int* n,
                const int* nb,
                int* kb,
                float* a,
                const int* lda,
                float* e,
                int* ipiv,
                float* w,
                const int* ldw,
                int* info);
            
void dlasyf_rk_(const char* uplo,
                const int* n,
                const int* nb,
                int* kb,
                double* a,
                const int* lda,
                double* e,
                int* ipiv,
                double* w,
                const int* ldw,
                int* info);

void clasyf_rk_(const char* uplo,
                const int* n,
                const int* nb,
                int* kb,
                kml_complex_float *a,
                const int* lda,
                kml_complex_float *e,
                int* ipiv,
                kml_complex_float *w,
                const int* ldw,
                int* info);

void zlasyf_rk_(const char* uplo,
                const int* n,
                const int* nb,
                int* kb,
                kml_complex_double *a,
                const int* lda,
                kml_complex_double *e,
                int* ipiv,
                kml_complex_double *w,
                const int* ldw,
                int* info);