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
    int threads_use = MIN(max_threads, 24);

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
        Xerbla(name, info, 8);
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
                LASYF_RK(uplo, &k, &nb, &kb, A, lda, E, ipiv, work, &ldwork,
                         &iInfo);
            } else {
                /*
                 *           Use unblocked code to factorize columns 1:k of A
                 */
                BlasSetNumThreads(MIN(max_threads, 4));
                SYTF2_RK(uplo, &k, A, lda, E, ipiv, &iInfo);
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
             *
             *        Apply permutations to the leading panel 1:k-1
             *
             *        Read IPIV from the last block factored, i.e.
             *        indices  k-kb+1:k and apply row permutations to the
             *        last k+1 colunms k+1:N after that block
             *        (We can do the simple loop over IPIV with decrement -1,
             *        since the ABS value of IPIV( I ) represents the row index
             *        of the interchange with row i in both 1x1 and 2x2 pivot
             * cases)
             */
            if (k < N) {
//#pragma omp parallel for private(ip, intTmp)
/*
for (i = k; i >= k - kb + 1; i--) {
    ip = ABS_(ipiv[i - 1]);
    if (ip != i) {
        intTmp = N - k;
        SWAP_(&intTmp, A + (i - 1) + (k) * (LDA), &LDA,
              A + (ip - 1) + (k) * (LDA), &LDA);
    }
}*/
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
            /*
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
                LASYF_RK(uplo, &param, &nb, &kb, A + k - 1 + (k - 1) * (LDA),
                         lda, E + k - 1, ipiv + k - 1, work, &ldwork, &iInfo);
            } else {
                /*
                 *           Use unblocked code to factorize columns k:n of A
                 */
                BlasSetNumThreads(MIN(max_threads, 4));
                int param = N - k + 1;
                SYTF2_RK(uplo, &param, A + k - 1 + (k - 1) * (LDA), lda,
                         E + k - 1, ipiv + k - 1, &iInfo);
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
             *        Apply permutations to the leading panel 1:k-1
             *
             *        Read IPIV from the last block factored, i.e.
             *        indices  k:k+kb-1 and apply row permutations to the
             *        first k-1 colunms 1:k-1 before that block
             *        (We can do the simple loop over IPIV with increment 1,
             *        since the ABS value of IPIV( I ) represents the row index
             *        of the interchange with row i in both 1x1 and 2x2 pivot
             * cases)
             */
            if (k > 1) {
                /*
                for (i = k; i <= k + kb - 1; i++) {
                    ip = ABS_(ipiv[i - 1]);
                    if (ip != i) {
                        intTmp = k - 1;
                        SWAP_(&intTmp, A + (i - 1), &LDA, A + (ip - 1), &LDA);
                    }
                }*/
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
    // End of SSYTRF_RK
}