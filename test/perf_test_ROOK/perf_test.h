#pragma once

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "../../include/SYSV_ROOK/SYSV_ROOK.h"
#include "../../include/SYSV_ROOK/type.h"
#include "creat_matrix.h"
#include "perT_type.h"

void perfTest(int N, int M, char uplo);

double calGflops(int N, int M, char uplo);

double _opla_(char* SUBNAM, int* M, int* N, int* KL, int* KU, int* NB);

void _sysv_rook_(const char* uplo,
                 const int* n,
                 const int* nrhs,
                 dataType* a,
                 const int* lda,
                 int* ipiv,
                 dataType* b,
                 const int* ldb,
                 dataType* work,
                 const int* lWork,
                 int* info);
