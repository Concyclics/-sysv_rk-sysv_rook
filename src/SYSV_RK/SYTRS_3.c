/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYTRS_3.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#include "SYTRS_3.h"

void SYTRS_3(const char* uplo,
             const int* n,
             const int* nrhs,
             dataType* A,
             const int* lda,
             dataType* E,
             int* ipiv,
             dataType* B,
             const int* ldb,
             int* info) {
    const dataType CONE = 1;
    const dataType NEG_CONE = -1;

    const int maxThreads = omp_get_max_threads();
    const int N = *n;
    const int NRHS = *nrhs;
    const int LDA = *lda;
    const int LDB = *ldb;
    int upper;
    int i, j, k, kp;
    dataType ak, akm1, akm1k, bk, bkm1, deNom;

    int swaped[N];
    int visited[N];
    int index[N];
    int cnt;

    *info = 0;
    upper = KmlLsame(*uplo, 'U');
    if (!upper && !KmlLsame(*uplo, 'L')) {
        *info = -1;
    } else if (N < 0) {
        *info = -2;
    } else if (NRHS < 0) {
        *info = -3;
    } else if (LDA < MAX(1, N)) {
        *info = -5;
    } else if (LDB < MAX(1, N)) {
        *info = -9;
    }
    int neg_info = -*info;
    if (*info != 0) {
#ifdef SINGLE
        Xerbla("SSYTRS_3", &neg_info, 8);
#endif
#ifdef DOUBLE
        Xerbla("DSYTRS_3", &neg_info, 8);
#endif
#ifdef COMPLEX
        Xerbla("CSYTRS_3", &neg_info, 8);
#endif
#ifdef COMPLEX16
        Xerbla("ZSYTRS_3", &neg_info, 8);
#endif

        return;
    }
    /*
     *     Quick return if possible
     */
    if (N == 0 || NRHS == 0) {
        return;
    }
    if (upper) {
        if (maxThreads == 1 && N < 1000) {
            for (k = N; k > 0; k--) {
                kp = ABS_(ipiv[k - 1]);
                if (kp != k) {
                    SWAP_(&NRHS, B + k - 1, &LDB, B + kp - 1, &LDB);
                }
            }
        } else {
            for (k = 0; k < N; k++) {
                swaped[k] = k;
                visited[k] = 0;
            }
            for (k = N - 1; k > 0; k--) {
                kp = ABS_(ipiv[k]);
                kp--;
                if (kp != k) {
                    int tmp = swaped[k];
                    swaped[k] = swaped[kp];
                    swaped[kp] = tmp;
                    visited[k] = 1;
                    visited[kp] = 1;
                }
            }

            cnt = 0;
            for (k = 0; k < N; k++) {
                if (visited[k] == 1) {
                    index[cnt] = k;
                    cnt++;
                }
            }

#pragma omp parallel for private(k, i)
            for (k = 0; k < NRHS; k++) {
                dataType Btmp[cnt];
                for (i = 0; i < cnt; i++) {
                    Btmp[i] = B[swaped[index[i]] + LDB * k];
                }
                for (i = 0; i < cnt; i++) {
                    B[index[i] + LDB * k] = Btmp[i];
                }
            }
        }
#if defined(COMPLEX) || defined(COMPLEX16)
        TRSM_("L", "U", "N", "U", &N, &NRHS, (void*)&CONE, A, &LDA, B, &LDB);
#else
        TRSM_("L", "U", "N", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
#endif

        int cnt_pos = 0;
        int cnt_neg = 0;
        int pos_num[N];
        int neg_num[N];
        for (i = N; i > 0; i--) {
            if (ipiv[i - 1] > 0) {
                pos_num[cnt_pos] = i;
                cnt_pos++;
            } else {
                neg_num[cnt_neg] = i;
                cnt_neg++;
                i--;
            }
        }

        int task_i;
#pragma omp parallel for private(i)
        for (task_i = 0; task_i < cnt_pos; task_i++) {
            i = pos_num[task_i];
            dataType ALPHA_ = T_ONE / A[i - 1 + (i - 1) * LDA];
#if defined(COMPLEX) || defined(COMPLEX16)
            SCAL_(&NRHS, (void*)&ALPHA_, B + i - 1, &LDB);
#else
            SCAL_(&NRHS, &ALPHA_, B + i - 1, &LDB);
#endif
        }

#pragma omp parallel for private(i, j, deNom, akm1k, akm1, ak, bkm1, bk)
        for (task_i = 0; task_i < cnt_neg; task_i++) {
            i = neg_num[task_i];
            akm1k = E[i - 1];
            akm1 = A[i - 2 + (i - 2) * LDA] / akm1k;
            ak = A[i - 1 + (i - 1) * LDA] / akm1k;
            deNom = akm1 * ak - T_ONE;
            for (j = 1; j <= NRHS; j++) {
                bkm1 = B[i - 2 + (j - 1) * LDB] / akm1k;
                bk = B[i - 1 + (j - 1) * LDB] / akm1k;
                B[i - 2 + (j - 1) * LDB] = (ak * bkm1 - bk) / deNom;
                B[i - 1 + (j - 1) * LDB] = (akm1 * bk - bkm1) / deNom;
            }
        }

#if defined(COMPLEX) || defined(COMPLEX16)
        TRSM_("L", "U", "T", "U", &N, &NRHS, (void*)&CONE, A, &LDA, B, &LDB);
#else
        TRSM_("L", "U", "T", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
#endif
        if (maxThreads == 1 && N < 1000) {
            for (k = 1; k <= N; k++) {
                kp = ABS_(ipiv[k - 1]);
                if (kp != k) {
                    SWAP_(&NRHS, B + k - 1, &LDB, B + kp - 1, &LDB);
                }
            }
        } else {
            for (k = 0; k < N; k++) {
                swaped[k] = k;
                visited[k] = 0;
            }

            for (k = 0; k < N; k++) {
                kp = ABS_(ipiv[k]);
                kp--;
                if (kp != k) {
                    int tmp = swaped[k];
                    swaped[k] = swaped[kp];
                    swaped[kp] = tmp;
                    visited[k] = 1;
                    visited[kp] = 1;
                }
            }

            cnt = 0;
            for (k = 0; k < N; k++) {
                if (visited[k] == 1) {
                    index[cnt] = k;
                    cnt++;
                }
            }

#pragma omp parallel for private(k, i)
            for (k = 0; k < NRHS; k++) {
                dataType Btmp[cnt];
                for (i = 0; i < cnt; i++) {
                    Btmp[i] = B[swaped[index[i]] + k * LDB];
                }
                for (i = 0; i < cnt; i++) {
                    B[index[i] + k * LDB] = Btmp[i];
                }
            }
        }
    } else {
        if (maxThreads == 1 && N < 1000) {
            for (k = 1; k <= N; k++) {
                kp = ABS_(ipiv[k - 1]);
                if (kp != k) {
                    SWAP_(&NRHS, B + k - 1, &LDB, B + kp - 1, &LDB);
                }
            }
        } else {
            for (k = 0; k < N; k++) {
                swaped[k] = k;
                visited[k] = 0;
            }

            for (k = 0; k < N; k++) {
                kp = ABS_(ipiv[k]);
                kp--;
                if (kp != k) {
                    int tmp = swaped[k];
                    swaped[k] = swaped[kp];
                    swaped[kp] = tmp;
                    visited[k] = 1;
                    visited[kp] = 1;
                }
            }

            cnt = 0;
            for (k = 0; k < N; k++) {
                if (visited[k] == 1) {
                    index[cnt] = k;
                    cnt++;
                }
            }

#pragma omp parallel for private(k, i)
            for (k = 0; k < NRHS; k++) {
                dataType Btmp[cnt];
                for (i = 0; i < cnt; i++) {
                    Btmp[i] = B[swaped[index[i]] + k * LDB];
                }
                for (i = 0; i < cnt; i++) {
                    B[index[i] + k * LDB] = Btmp[i];
                }
            }
        }
#if defined(COMPLEX) || defined(COMPLEX16)
        TRSM_("L", "L", "N", "U", &N, &NRHS, (void*)&CONE, A, &LDA, B, &LDB);
#else
        TRSM_("L", "L", "N", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
#endif

        int cnt_pos = 0;
        int cnt_neg = 0;
        int pos_num[N];
        int neg_num[N];
        for (i = 1; i <= N; i++) {
            if (ipiv[i - 1] > 0) {
                pos_num[cnt_pos] = i;
                cnt_pos++;
            } else {
                neg_num[cnt_neg] = i;
                cnt_neg++;
                i++;
            }
        }
        int task_i;
#pragma omp parallel for private(i)
        for (task_i = 0; task_i < cnt_pos; task_i++) {
            i = pos_num[task_i];
            dataType ALPHA_ = T_ONE / A[i - 1 + (i - 1) * LDA];
#if defined(COMPLEX) || defined(COMPLEX16)
            SCAL_(&NRHS, (void*)&ALPHA_, B + i - 1, &LDB);
#else
            SCAL_(&NRHS, &ALPHA_, B + i - 1, &LDB);
#endif
        }

#pragma omp parallel for private(i, j, deNom, akm1k, akm1, ak, bkm1, bk)
        for (task_i = 0; task_i < cnt_neg; task_i++) {
            i = neg_num[task_i];
            akm1k = E[i - 1];
            akm1 = A[i - 1 + (i - 1) * LDA] / akm1k;
            ak = A[i + i * LDA] / akm1k;
            deNom = akm1 * ak - T_ONE;
            for (j = 1; j <= NRHS; j++) {
                bkm1 = B[i - 1 + (j - 1) * LDB] / akm1k;
                bk = B[i + (j - 1) * LDB] / akm1k;
                B[i - 1 + (j - 1) * LDB] = (ak * bkm1 - bk) / deNom;
                B[i + (j - 1) * LDB] = (akm1 * bk - bkm1) / deNom;
            }
        }
#if defined(COMPLEX) || defined(COMPLEX16)
        TRSM_("L", "L", "T", "U", &N, &NRHS, (void*)&CONE, A, &LDA, B, &LDB);
#else
        TRSM_("L", "L", "T", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
#endif
        if (maxThreads == 1 && N < 1000) {
            for (k = N; k >= 1; k--) {
                kp = ABS_(ipiv[k - 1]);
                if (kp != k) {
                    SWAP_(&NRHS, B + k - 1, &LDB, B + kp - 1, &LDB);
                }
            }
        } else {
            for (k = 0; k < N; k++) {
                swaped[k] = k;
                visited[k] = 0;
            }

            for (k = N - 1; k >= 0; k--) {
                kp = ABS_(ipiv[k]);
                kp--;
                if (kp != k) {
                    int tmp = swaped[k];
                    swaped[k] = swaped[kp];
                    swaped[kp] = tmp;
                    visited[k] = 1;
                    visited[kp] = 1;
                }
            }

            cnt = 0;
            for (k = 0; k < N; k++) {
                if (visited[k] == 1) {
                    index[cnt] = k;
                    cnt++;
                }
            }

#pragma omp parallel for private(k, i)
            for (k = 0; k < NRHS; k++) {
                dataType Btmp[cnt];
                for (i = 0; i < cnt; i++) {
                    Btmp[i] = B[swaped[index[i]] + k * LDB];
                }
                for (i = 0; i < cnt; i++) {
                    B[index[i] + k * LDB] = Btmp[i];
                }
            }
        }
    }
    return;
}