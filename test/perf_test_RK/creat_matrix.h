#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "funT_type.h"

/**
 * ////////////////////////////////abs///////////////////////////////////////
 * @Creat a random symmetry matrix A in size N*N and a random matrix B in size
 * N*NRHS, and the data type is float.
 * @param[in]		uplo		uplo='L'or'U' indicates A is stored as a
 * lower triangular matrix or a upper triangular matrix.
 * @param[in]		N
 * @param[in]		NRHS
 * @param[in out]	A			A random symmetry matrix in size
 * N*N.
 * @param[in out]	B			A random symmetry matrix in size
 * N*NRHS.
 * @retva
 */
void _CreatMatrix(char uplo,
                  int N,
                  int NRHS,
                  funcTdataType* A,
                  funcTdataType* B);
