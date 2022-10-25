/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYSV_RK.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#include "SYSV_RK.h"

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
             int* info) {
    int lwkopt;
    /*
     *     Test the input parameters.
     */
    if (info == NULL) {
        int neg_info = 12;
#ifdef SINGLE
        Xerbla("ssysv_rk", &neg_info, 12);
#endif
#ifdef DOUBLE
        Xerbla("dsysv_rk", &neg_info, 12);
#endif
#ifdef COMPLEX
        Xerbla("csysv_rk", &neg_info, 12);
#endif
#ifdef COMPLEX16
        Xerbla("zsysv_rk", &neg_info, 12);
#endif
        return;
    }
    *info = 0;
    if (uplo == NULL || ((!KmlLsame(*uplo, 'U')) && (!KmlLsame(*uplo, 'L')))) {
        *info = -1;
    } else if (n == NULL || *n < 0) {
        *info = -2;
    } else if (nrhs == NULL || *nrhs < 0) {
        *info = -3;
    } else if (a == NULL) {
        *info = -4;
    } else if (lda == NULL || *lda < MAX(1, *n)) {
        *info = -5;
    } else if (e == NULL) {
        *info = -6;
    } else if (ipiv == NULL) {
        *info = -7;
    } else if (b == NULL) {
        *info = -8;
    } else if (ldb == NULL || *ldb < MAX(1, *n)) {
        *info = -9;
    } else if (work == NULL) {
        *info = -10;
    } else if (lwork == NULL || (*lwork < 1 && *lwork != -1)) {
        *info = -11;
    }

    if (*info == 0) {
        if (*n != 0) {
            lwkopt = 1;
        } else {
            int temp = -1;
            SYTRF_RK(uplo, n, a, lda, e, ipiv, work, &temp, info);
            lwkopt = work[0];
        }
        work[0] = lwkopt;
    }

    if (*info != 0) {
        int neg_info = -*info;
#ifdef SINGLE
        Xerbla("ssysv_rk", &neg_info, 12);
#endif
#ifdef DOUBLE
        Xerbla("dsysv_rk", &neg_info, 12);
#endif
#ifdef COMPLEX
        Xerbla("csysv_rk", &neg_info, 12);
#endif
#ifdef COMPLEX16
        Xerbla("zsysv_rk", &neg_info, 12);
#endif
        return;
    } else if (*lwork == -1) {
        return;
    }

    /*
     *     Compute the factorization A = P*U*D*(U**T)*(P**T) or
     *     A = P*U*D*(U**T)*(P**T).
     */

    SYTRF_RK(uplo, n, a, lda, e, ipiv, work, lwork, info);
    if (*info == 0) {
        // Solve the system A*X = B with BLAS3 solver, overwriting B with X.
        SYTRS_3(uplo, n, nrhs, a, lda, e, ipiv, b, ldb, info);
    }

    work[0] = lwkopt;
    return;
}