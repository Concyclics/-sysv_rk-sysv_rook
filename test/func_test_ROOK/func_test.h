/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of _FuncTest, _pot02_ and _sysv_rook_.
 * Author: linhouzhong
 * Create: 2022-06-26
 */

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include "../../include/SYSV_ROOK/SYSV_ROOK.h"
#include "creat_matrix.h"
#include "funT_type.h"
/**
 * ////////////////////////////////_FuncTest///////////////////////////////////////
 * @Brief Fuction test.
 * @param[in]		timess		Times of result test.
 * @param[in]		scale		scale of result test.
 * @retva
 */
void _FuncTest(int times, char scale);

int cover();

double _TestWithSeed(int rand_seed,
                     int N,
                     int NRHS,
                     char uplo,
                     funcTdataType* A,
                     funcTdataType* A1,
                     funcTdataType* B1,
                     funcTdataType* B3,
                     funcTdataType* W,
                     funcTdataAccu* W1,
                     int* ipiv);
/**
 * ////////////////////////////////_pot02_///////////////////////////////////////
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
 * @param[in]	    a		    The original symmetric matrix A.
 * @param[in]		lda		    The leading dimension of the array
 A.
 * @param[in]		x		    The computed solution vectors for
 the system of linear equations.
 * @param[in]		ldx		    The leading dimension of the array
 X.
 * @param[in,out]	b		    On entry, the right hand side
 vectors for the system of linear equations. On exit, B is overwritten with the
 difference B - A*X.
 * @param[in]		ldb		    The leading dimension of the array
 B.
 * @param[in]		rwork		The residual work array.
 * @param[out]		resid		The maximum over the number of right
 hand sides of norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
 * @retva
 */
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

/**
 * ////////////////////////////////_sysv_rook_///////////////////////////////////////
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
void _sysv_rook_(const char* uplo,
                 const int* n,
                 const int* nrhs,
                 funcTdataType* a,
                 const int* lda,
                 int* ipiv,
                 funcTdataType* b,
                 const int* ldb,
                 funcTdataType* work,
                 const int* lWork,
                 int* info);

/**
 * ////////////////////////////////nullptr_test///////////////////////////////////////
 * @Brief null pointer test.
 * @param[in]		N		A is a N*N matrix.
 * @param[in]		M		B is a N*M matrix.
 * @retva
 */
void nullptr_test(int N, int M);

void exception_test(int N, int M);