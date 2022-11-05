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

void LASYF_ROOK(const char* uplo,
                const int* n,
                const int* nb,
                int* kb,
                dataType* a,
                const int* lda,
                int* ipiv,
                dataType* w,
                const int* ldw,
                int* info);
