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

void SYTRF_RK(const char* uplo,
              const int* n,
              dataType* A,
              const int* lda,
              dataType* E,
              int* ipiv,
              dataType* work,
              const int* lwork,
              int* info);