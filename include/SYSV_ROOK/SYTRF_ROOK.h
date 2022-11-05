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

void SYTRF_ROOK(const char* uplo,
                const int* n,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* work,
                const int* lwork,
                int* info);