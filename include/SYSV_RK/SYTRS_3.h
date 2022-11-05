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