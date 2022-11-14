/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of _FuncTest, _pot02_ and _sysv_rk_.
 * Author: linhouzhong
 * Create: 2022-06-26
 */

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "../../include/SYSV_RK/SYSV_RK.h"
#include "creat_matrix.h"
#include "funT_type.h"
/**
 * ////////////////////////////////_FuncTest///////////////////////////////////////
 * @Brief Fuction test.
 * @param[in]		timess		Times of result test.
 * @param[in]		scale		scale of result test.
 * @retva
 */
void _FuncTest(int times,
               char scale,
               funcTdataType* A,
               funcTdataType* A1,
               funcTdataType* E,
               funcTdataType* B1,
               funcTdataType* B2,
               funcTdataType* B3,
               funcTdataType* B4,
               funcTdataType* W,
               funcTdataAccu* W1,
               int* ipiv);

int cover();

double _TestWithSeed(int rand_seed,
                     int N,
                     int NRHS,
                     char uplo,
                     funcTdataType* A,
                     funcTdataType* A1,
                     funcTdataType* E,
                     funcTdataType* B1,
                     funcTdataType* B3,
                     funcTdataType* W,
                     funcTdataAccu* W1,
                     int* ipiv);

void _pot02_(const char* uplo,
             const int* n,
             const int* nrhs,
             funcTdataType* a,
             const int* lda,
             funcTdataType* x,
             const int* ldx,
             funcTdataType* b,
             const int* ldb,
             funcTdataAccu* rwork,
             funcTdataAccu* resid);

void _sysv_rk_(const char* uplo,
               const int* n,
               const int* nrhs,
               funcTdataType* a,
               const int* lda,
               funcTdataType* e,
               int* ipiv,
               funcTdataType* b,
               const int* ldb,
               funcTdataType* work,
               const int* lWork,
               int* info);

void nullptr_test(int N, int M);

void exception_test(int N, int M);
