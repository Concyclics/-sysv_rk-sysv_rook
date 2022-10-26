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

/**
 * ////////////////////////////////SYSV_RK///////////////////////////////////////
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
 * @param[out]		E		    Contains the output computed by the
 factorization routine DSYTRF_RK, i.e. the superdiagonal (or subdiagonal)
                                elements of the symmetric block diagonal matrix
 D with 1-by-1 or 2-by-2 diagonal blocks
 * @param[out]		ipiv		Details of the interchanges and the
 block structure of D, as determined by ?SYTRF_RK.
 * @param[in,out]	b		    On entry, the N-by-NRHS right hand
 side matrix B. On exit, if INFO = 0, the N-by-NRHS solution matrix X.
 * @param[in]		ldb		    The leading dimension of the array
 B.
 * @param[out]		work		Work array used in the factorization
 stage.
 * @param[in]		lWork		The length of WORK.
 * @param[out]		info		= 0: successful exit
                                < 0: If INFO = -k, the k-th argument had an
 illegal value > 0: If INFO = k, the matrix A is singular,
 */
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