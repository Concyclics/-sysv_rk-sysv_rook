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

/**
 * ////////////////////////////////LASYF_RK///////////////////////////////////////
 * @Brief Computes a partial factorization of a real symmetric matrix A using
 the bounded Bunch-Kaufman (rook) diagonal pivoting method.
 * @param[in]		uplo	    Specifies whether the upper or lower
 triangular part of the symmetric matrix A is stored: = 'U':  Upper triangle of
 A is stored; = 'L':  Lower triangle of A is stored.
 * @param[in]		n		    The number of linear equations,
 i.e., the order of the matrix A.
 * @param[in]		nb	        The maximum number of columns of the
 matrix A that should be factored.
 * @param[in]		kb	        The number of columns of A that were
 actually factored.
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
 * @param[out]		w		    Work array used in the factorization
 stage.
 * @param[in]		ldw		    The length of WORK.
 * @param[out]		info		= 0: successful exit
                                < 0: If INFO = -k, the k-th argument had an
 illegal value > 0: If INFO = k, the matrix A is singular,
 */
void LASYF_RK(const char* uplo,
              const int* n,
              const int* nb,
              int* kb,
              dataType* a,
              const int* lda,
              dataType* e,
              int* ipiv,
              dataType* w,
              const int* ldw,
              int* info);
