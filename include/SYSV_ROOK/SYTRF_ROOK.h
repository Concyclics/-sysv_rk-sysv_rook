/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of SYTRF_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-29
 *******************************************************************************/

#pragma once

#include <complex.h>
#include <math.h>
#include "../kml/kblas.h"
#include "../kml/kml_service.h"
#include "LASYF_ROOK.h"
#include "SYTF2_ROOK.h"
#include "ilaenv.h"
#include "type.h"

/**
 * ////////////////////////////////SYTRF_ROOK///////////////////////////////////////
 * @Brief Computes the residual for the solution of a symmetric system of linear
 equations  A*x = b: RESID = norm(B - A*X) / ( norm(A) * norm(X) * EPS ), where
 EPS is the machine epsilon.
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
 * @param[out]		work		Work array used in the factorization
 stage.
 * @param[in]		lWork		The length of WORK.
 * @param[out]		info		= 0: successful exit
                                < 0: If INFO = -k, the k-th argument had an
 illegal value > 0: If INFO = k, the matrix A is singular,
 */

void SYTRF_ROOK(const char* uplo,
                const int* n,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* work,
                const int* lwork,
                int* info);