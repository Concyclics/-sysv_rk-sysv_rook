/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYSV_ROOK.
 * Author: linhouzhong
 * Create: 2022-06-30
 *******************************************************************************/

#pragma once
#include <kml_service.h>
#include "SYTRF_ROOK.h"
#include "SYTRS_ROOK.h"
#include "type.h"

void SYSV_ROOK(const char* uplo,
               const int* n,
               const int* nrhs,
               dataType* a,
               const int* lda,
               int* ipiv,
               dataType* b,
               const int* ldb,
               dataType* work,
               const int* lwork,
               int* info);