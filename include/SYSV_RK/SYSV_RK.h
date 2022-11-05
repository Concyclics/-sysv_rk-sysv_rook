/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYSV_RK.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#pragma once

#include <kml_service.h>
#include "SYTRF_RK.h"
#include "SYTRS_3.h"
#include "type.h"

void SYSV_RK(const char* uplo,
             const int* n,
             const int* nrhs,
             dataType* a,
             const int* lda,
             dataType* e,
             int* ipiv,
             dataType* b,
             const int* ldb,
             dataType* work,
             const int* lwork,
             int* info);