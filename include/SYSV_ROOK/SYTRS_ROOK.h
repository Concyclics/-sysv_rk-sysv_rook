/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTRS_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-27
 *******************************************************************************/

#pragma once
#include <complex.h>
#include <math.h>
#include "../kml/kblas.h"
#include "../kml/kml_service.h"
#include "type.h"

/**
 * ////////////////////////////////SYTRS_ROOK//////////////////////////////////////
 * @Brief SSYTRS_ROOK solves a system of linear equations A * X = B with a
 symmetric matrix A using the factorization computed by SYTRF_ROOK: RESID =
 norm(B - A*X) / ( norm(A) * norm(X) * EPS ), where EPS is the machine epsilon.
 * @param[in]		uplo	    Specifies whether the upper or lower
 triangular part of the symmetric matrix A is stored: = 'U':  Upper triangle of
 A is stored; = 'L':  Lower triangle of A is stored.
 * @param[in]		n		    The number of linear equations,
 i.e., the order of the matrix A.
 * @param[in]		nrhs	    The number of right hand sides, i.e., the
 number of columns of the matrix B.
 * @param[in,out]	a		    On entry,the original symmetric
 matrix A.
 *                              On exit, if INFO = 0, diagonal of the block
 diagonal matrix D and factors U or L as computed by ?SYTRF_RK
 * @param[in]		lda		    The leading dimension of the array
 A.
 * @param[out]		ipiv		Details of the interchanges and the
 block structure of D, as determined by ?SYTRF_RK.
 * @param[in,out]	b		    On entry, the N-by-NRHS right hand
 side matrix B. On exit, if INFO = 0, the N-by-NRHS solution matrix X.
 * @param[in]		ldb		    The leading dimension of the array
 B.
 * @param[out]		info		= 0: successful exit
                                < 0: If INFO = -k, the k-th argument had an
 illegal value > 0: If INFO = k, the matrix A is singular,
 */

void SYTRS_ROOK(const char* uplo,
                const int* n,
                const int* nrhs,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* B,
                const int* ldb,
                int* info);