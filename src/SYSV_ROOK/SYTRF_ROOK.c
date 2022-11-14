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
    int threads_use;

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
        int neg_info = -*info;
#ifdef SINGLE
        Xerbla("SYTRF_ROOK", &neg_info, 8);
#endif
#ifdef DOUBLE
        Xerbla("DYTRF_ROOK", &neg_info, 8);
#endif
#ifdef COMPLEX
        Xerbla("CYTRF_ROOK", &neg_info, 8);
#endif
#ifdef COMPLEX16
        Xerbla("ZYTRF_ROOK", &neg_info, 8);
#endif
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
        k = N;
        while (k > 0) {
            if (k > nb) {
                if (N <= 8500) {
                    threads_use = 16;
                } else if (N <= 16500) {
                    threads_use = 20;
                } else if (N <= 25000) {
                    threads_use = 24;
                } else if (N <= 30000) {
                    threads_use = 28;
                } else if (N <= 40000) {
                    threads_use = 32;
                } else {
                    threads_use = 48;
                }
                BlasSetNumThreads(MIN(max_threads, threads_use));
                LASYF_ROOK(uplo, &k, &nb, &kb, A, lda, ipiv, work, &ldWork,
                           &iInfo);
            } else {
                BlasSetNumThreads(MIN(max_threads, 4));
                SYTF2_ROOK(uplo, &k, A, lda, ipiv, &iInfo);
                kb = k;
            }
            if (info == 0 && iInfo > 0) {
                *info = iInfo;
            }
            k -= kb;
        }
    } else {
        k = 1;
        while (k <= N) {
            if (k <= N - nb) {
                if (N <= 8500) {
                    threads_use = 16;
                } else if (N <= 16500) {
                    threads_use = 20;
                } else if (N <= 25000) {
                    threads_use = 24;
                } else if (N <= 30000) {
                    threads_use = 28;
                } else if (N <= 40000) {
                    threads_use = 32;
                } else {
                    threads_use = 48;
                }
                BlasSetNumThreads(MIN(max_threads, threads_use));
                int param = N - k + 1;
                LASYF_ROOK(uplo, &param, &nb, &kb, A + k - 1 + (k - 1) * (LDA),
                           lda, ipiv + k - 1, work, &ldWork, &iInfo);
            } else {
                BlasSetNumThreads(MIN(max_threads, 4));
                int param = N - k + 1;
                SYTF2_ROOK(uplo, &param, A + k - 1 + (k - 1) * (LDA), lda,
                           ipiv + k - 1, &iInfo);
                kb = N - k + 1;
            }
            if (info == 0 && iInfo > 0) {
                *info = iInfo + k - 1;
            }
            for (i = k; i < k + kb; i++) {
                if (ipiv[i - 1] > 0) {
                    ipiv[i - 1] += k - 1;
                } else {
                    ipiv[i - 1] += -k + 1;
                }
            }
            k += kb;
        }
    }
    work[0] = (float)lwkopt;
    BlasSetNumThreads(max_threads);
    return;
}