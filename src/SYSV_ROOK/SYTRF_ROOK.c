/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYTRF_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-29
 *******************************************************************************/

#include "SYTRF_ROOK.h"

void SYTRF_ROOK(const char* uplo,
                const int* n,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* work,
                const int* lwork,
                int* info) {
    const int N = *n;
    const int LDA = *lda;
    const int LWORK = *lwork;
    int intTmp;

    int lQuery, upper;
    int i, iInfo, ip, iws, k, kb, ldWork, lwkopt, nb, nbMin;
    int one1 = 1;
    int _one_1 = -1;
    int two2 = 2;

    int max_threads = omp_get_max_threads();
    int threads_use = MIN(max_threads, 24);

    *info = 0;
    upper = KmlLsame(*uplo, 'U');
    lQuery = (*lwork == -1);
    if (!upper && !KmlLsame(*uplo, 'L')) {
        *info = -1;
    } else if (N < 0) {
        *info = -2;
    } else if (LDA < MAX(1, N)) {
        *info = -4;
    } else if (LWORK < 1 && !lQuery) {
        *info = -7;
    }

    if (*info != 0) {
        Xerbla("SYTRF_ROOK", info, 8);
        return;
    } else if (lQuery) {
        return;
    }

    if (*info == 0) {
        nb = ilaenv(&one1, "SYTRF_ROOK", uplo, n, &_one_1, &_one_1, &_one_1);
        lwkopt = MAX(1, (N)*nb);
        work[0] = (float)lwkopt;
    }

    if (N == 0) {
        work[0] = 1.0;
        return;
    }

    nbMin = 2;
    ldWork = N;
    if (nb > 1 && nb < N) {
        iws = ldWork * nb;
        if (LWORK < iws) {
            nb = MAX(LWORK / ldWork, 1);
            nbMin = MAX(2, ilaenv(&two2, "SYTRF_ROOK", uplo, n, &_one_1,
                                  &_one_1, &_one_1));
        }
    } else {
        iws = 1;
    }

    if (nb < nbMin) {
        nb = N;
    }

    if (upper) {
        /*
         *        Factorize A as U*D*U**T using the upper triangle of A
         *
         *        K is the main loop index, decreasing from N to 1 in steps of
         *        KB, where KB is the number of columns factorized by SLASYF_RK;
         *        KB is either NB or NB-1, or K for the last block
         */
        k = N;
        /*
         *        If K < 1, exit from loop
         */
        while (k > 0) {
            if (k > nb) {
                /*
                 *           Factorize columns k-kb+1:k of A and use blocked
                 * code to update columns 1:k-kb
                 */
                BlasSetNumThreads(threads_use);
                LASYF_ROOK(uplo, &k, &nb, &kb, A, lda, ipiv, work, &ldWork,
                           &iInfo);
            } else {
                /*
                 *           Use unblocked code to factorize columns 1:k of A
                 */
                BlasSetNumThreads(MIN(max_threads, 4));
                SYTF2_ROOK(uplo, &k, A, lda, ipiv, &iInfo);
                kb = k;
            }
            /*
             *        Set INFO on the first occurrence of a zero pivot
             */
            if (info == 0 && iInfo > 0) {
                *info = iInfo;
            }
            /*
             *        No need to adjust IPIV
             *
             *        Decrease K and return to the start of the main loop
             */
            k -= kb;
        }
        /*
         *        This label is the exit from main loop over K decreasing
         *        from N to 1 in steps of KB
         */
    } else {
        /*
         *        Factorize A as L*D*L**T using the lower triangle of A
         *
         *        K is the main loop index, increasing from 1 to N in steps of
         *        KB, where KB is the number of columns factorized by SLASYF_RK;
         *        KB is either NB or NB-1, or N-K+1 for the last block
         */
        k = 1;
        while (k <= N) {
            if (k <= N - nb) {
                /*
                 *           Factorize columns k:k+kb-1 of A and use blocked
                 * code to update columns k+kb:n
                 *
                 */
                BlasSetNumThreads(threads_use);
                int param = N - k + 1;
                LASYF_ROOK(uplo, &param, &nb, &kb, A + k - 1 + (k - 1) * (LDA),
                           lda, ipiv + k - 1, work, &ldWork, &iInfo);
            } else {
                /*
                 *           Use unblocked code to factorize columns k:n of A
                 */
                BlasSetNumThreads(MIN(max_threads, 4));
                int param = N - k + 1;
                SYTF2_ROOK(uplo, &param, A + k - 1 + (k - 1) * (LDA), lda,
                           ipiv + k - 1, &iInfo);
                kb = N - k + 1;
            }
            /*
             *        Set INFO on the first occurrence of a zero pivot
             */
            if (info == 0 && iInfo > 0) {
                *info = iInfo + k - 1;
            }
            /*
             *        Adjust IPIV
             */
            for (i = k; i < k + kb; i++) {
                if (ipiv[i - 1] > 0) {
                    ipiv[i - 1] += k - 1;
                } else {
                    ipiv[i - 1] += -k + 1;
                }
            }
            /*
             *        Increase K and return to the start of the main loop
             */
            k += kb;
        }
        /*
         *        This label is the exit from main loop over K increasing
         *        from 1 to N in steps of KB
         */

        /*
         *     End Lower
         */
    }
    work[0] = (float)lwkopt;
    BlasSetNumThreads(max_threads);
    return;
    // End of SSYTRF_ROOK
}