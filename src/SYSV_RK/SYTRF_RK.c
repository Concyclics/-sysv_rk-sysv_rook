/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYTRF_RK.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#include "SYTRF_RK.h"

void SYTRF_RK(const char* uplo,
              const int* n,
              dataType* A,
              const int* lda,
              dataType* E,
              int* ipiv,
              dataType* work,
              const int* lwork,
              int* info) {
    const int N = *n;
    const int LDA = *lda;
    const int LWORK = *lwork;
    int intTmp;
    int swaped[N];
    int visited[N];

    int lQuery, upper;
    int i, iInfo, ip, iws, k, kb, ldwork, lwkopt, nb, nbmin;
    int one1 = 1;
    int _one_1 = -1;
    int two2 = 2;

    int max_threads = omp_get_max_threads();
    // int threads_use = MIN(max_threads, 24);
    int threads_use;

    char* name;
    char* opts;
#ifdef SINGLE
    name = "SSYTRF_RK";
#elif defined DOUBLE
    name = "DSYTRF_RK";
#elif defined COMPLEX
    name = "CSYTRF_RK";
#elif defined COMPLEX16
    name = "ZSYTRF_RK";
#endif

    *info = 0;
    upper = KmlLsame(*uplo, 'U');
    lQuery = (LWORK == -1);
    if (!upper && !KmlLsame(*uplo, 'L')) {
        *info = -1;
    } else if (N < 0) {
        *info = -2;
    } else if (LDA < MAX(1, N)) {
        *info = -4;
    } else if (LWORK < 1 && !lQuery) {
        *info = -8;
    }

    if (*info != 0) {
        int negInfo = -*info;
        Xerbla(name, &negInfo, 8);
        return;
    } else if (lQuery) {
        return;
    }

    if (*info == 0) {
        nb = ilaenv(&one1, name, "RK", n, &_one_1, &_one_1, &_one_1);
        lwkopt = (N)*nb;
        work[0] = (float)lwkopt;
    }

    if (N == 0) {
        work[0] = 1.0;
        return;
    }

    nbmin = 2;
    ldwork = N;
    if (nb > 1 && nb < N) {
        iws = ldwork * nb;
        if (LWORK < iws) {
            nb = MAX(LWORK / ldwork, 1);
            nbmin =
                MAX(2, ilaenv(&two2, name, "RK", n, &_one_1, &_one_1, &_one_1));
        }
    } else {
        iws = 1;
    }

    if (nb < nbmin) {
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
                LASYF_RK(uplo, &k, &nb, &kb, A, lda, E, ipiv, work, &ldwork,
                         &iInfo);
            } else {
                BlasSetNumThreads(MIN(max_threads, 4));
                SYTF2_RK(uplo, &k, A, lda, E, ipiv, &iInfo);
                kb = k;
            }
            if (info == 0 && iInfo > 0) {
                *info = iInfo;
            }
            if (k < N) {
#pragma omp parallel for
                for (i = 0; i < N; i++) {
                    swaped[i] = i;
                    visited[i] = 0;
                }
                for (i = k - 1; i >= k - kb; i--) {
                    ip = ABS_(ipiv[i]);
                    ip--;
                    if (ip != i) {
                        int tmp;
                        tmp = swaped[ip];
                        swaped[ip] = swaped[i];
                        swaped[i] = tmp;
                        visited[ip] = 1;
                        visited[i] = 1;
                    }
                }
                int cnt = 0;
                int index[kb << 1];
                for (i = 0; i < N; i++) {
                    if (visited[i] == 1) {
                        index[cnt] = i;
                        cnt++;
                    }
                }

#pragma omp parallel for private(i)
                for (i = k; i < N; i++) {
                    dataType Atmp[cnt];
                    int j;
                    for (j = 0; j < cnt; j++) {
                        Atmp[j] = A[swaped[index[j]] + i * LDA];
                    }
                    for (j = 0; j < cnt; j++) {
                        A[index[j] + i * LDA] = Atmp[j];
                    }
                }
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
                LASYF_RK(uplo, &param, &nb, &kb, A + k - 1 + (k - 1) * (LDA),
                         lda, E + k - 1, ipiv + k - 1, work, &ldwork, &iInfo);
            } else {
                BlasSetNumThreads(MIN(max_threads, 4));
                int param = N - k + 1;
                SYTF2_RK(uplo, &param, A + k - 1 + (k - 1) * (LDA), lda,
                         E + k - 1, ipiv + k - 1, &iInfo);
                kb = N - k + 1;
            }
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
            if (k > 1) {
                for (i = 0; i < N; i++) {
                    swaped[i] = i;
                    visited[i] = 0;
                }
                for (i = k - 1; i < k + kb - 1; i++) {
                    ip = ABS_(ipiv[i]);
                    ip--;
                    if (ip != i) {
                        int tmp;
                        tmp = swaped[ip];
                        swaped[ip] = swaped[i];
                        swaped[i] = tmp;
                        visited[ip] = 1;
                        visited[i] = 1;
                    }
                }
                int cnt = 0;
                int index[kb << 1];
                for (i = 0; i < N; i++) {
                    if (visited[i] == 1) {
                        index[cnt] = i;
                        cnt++;
                    }
                }

#pragma omp parallel for private(i)
                for (i = 0; i < k - 1; i++) {
                    dataType Atmp[cnt];
                    int j;
                    for (j = 0; j < cnt; j++) {
                        Atmp[j] = A[swaped[index[j]] + i * LDA];
                    }
                    for (j = 0; j < cnt; j++) {
                        A[index[j] + i * LDA] = Atmp[j];
                    }
                }
            }
            k += kb;
        }
    }
    work[0] = (float)lwkopt;
    BlasSetNumThreads(max_threads);
    return;
}