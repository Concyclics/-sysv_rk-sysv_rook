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

void SYTRS_ROOK(const char* uplo,
                const int* n,
                const int* nrhs,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* B,
                const int* ldb,
                int* info);