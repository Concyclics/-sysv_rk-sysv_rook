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
        Xerbla("SSYTRS_3", &neg_info, 10);
    }
    /*
     *     Quick return if possible
     */
    if (N == 0 || NRHS == 0) {
        return;
    }
    if (upper) {
/*
 *        Begin Upper
 *
 *        Solve A*X = B, where A = U*D*U**T.
 *
 *        P**T * B
 *
 *        Interchange rows K and IPIV(K) of matrix B in the same order
 *        that the formation order of IPIVi vector for Upper case.
 *
 *        (We can do the simple loop over IPIV with decrement -1,
 *        since the ABS value of IPIVi represents the row index
 *        of the interchange with row i in both 1x1 and 2x2 pivot cases)
 */
//#pragma omp parallel for schedule(static) private(k, kp)
/*
//#pragma omp parallel for private(k, kp)
for (k = N; k > 0; k--) {
    kp = ABS_(ipiv[k - 1]);
    if (kp != k) {
        //printf("k = %d, kp = %d\n", k, kp);
        SWAP_(&NRHS, B + k - 1, &LDB, B + kp - 1, &LDB);
    }
}
*/
// new swap
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
        for (k = 0; k < LDB; k++) {
            dataType Btmp[cnt];
            for (i = 0; i < cnt; i++) {
                Btmp[i] = B[swaped[index[i]] + LDB * k];
            }
            for (i = 0; i < cnt; i++) {
                B[index[i] + LDB * k] = Btmp[i];
            }
        }

        /*
         *        Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
         */

        TRSM_("L", "U", "N", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
        /*
         *        Compute D \ B -> B   [ D \ (U \P**T * B) ]
         */

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

        /*
        for (i = N; i > 0; i--) {
            if (ipiv[i - 1] > 0) {
                dataType ALPHA_ = T_ONE / A[i - 1 + (i - 1) * LDA];
            #if defined(COMPLEX) || defined(COMPLEX16)
                SCAL_(&NRHS, (void*)&ALPHA_, B + i - 1, &LDB);
            #else
                SCAL_(&NRHS, &ALPHA_, B + i - 1, &LDB);
            #endif
            } else if (i > 1) {
                akm1k = E[i - 1];
                akm1 = A[i - 2 + (i - 2) * LDA] / akm1k;
                ak = A[i - 1 + (i - 1) * LDA] / akm1k;
                deNom = akm1 * ak - T_ONE;
            #pragma omp parallel for schedule(static) private(j)
                for (j = 1; j <= NRHS; j++) {
                    bkm1 = B[i - 2 + (j - 1) * LDB] / akm1k;
                    bk = B[i - 1 + (j - 1) * LDB] / akm1k;
                    B[i - 2 + (j - 1) * LDB] = (ak * bkm1 - bk) / deNom;
                    B[i - 1 + (j - 1) * LDB] = (akm1 * bk - bkm1) / deNom;
                }
                i--;
            }
        }*/

        /*
         *        Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
         */
        TRSM_("L", "U", "T", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
/*
 *        P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
 *
 *        Interchange rows K and IPIV(K) of matrix B in reverse order
 *        from the formation order of IPIVi vector for Upper case.
 *
 *        (We can do the simple loop over IPIV with increment 1,
 *        since the ABS value of IPIVi represents the row index
 *        of the interchange with row i in both 1x1 and 2x2 pivot cases)
 */
//#pragma omp parallel for schedule(static) private(k, kp)
//#pragma omp parallel for private(k, kp)
/*
for (k = 1; k <= N; k++) {
    kp = ABS_(ipiv[k - 1]);
    if (kp != k) {
        SWAP_(&NRHS, B + k - 1, &LDB, B + kp - 1, &LDB);
    }
}*/
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
        for (k = 0; k < LDB; k++) {
            dataType Btmp[cnt];
            for (i = 0; i < cnt; i++) {
                Btmp[i] = B[swaped[index[i]] + k * LDB];
            }
            for (i = 0; i < cnt; i++) {
                B[index[i] + k * LDB] = Btmp[i];
            }
        }

    } else {
/*
 *        Begin Lower
 *
 *        Solve A*X = B, where A = L*D*L**T.
 *
 *        P**T * B
 *        Interchange rows K and IPIV(K) of matrix B in the same order
 *        that the formation order of IPIVi vector for Lower case.
 *
 *        (We can do the simple loop over IPIV with increment 1,
 *        since the ABS value of IPIVi represents the row index
 *        of the interchange with row i in both 1x1 and 2x2 pivot cases)
 */
//#pragma omp parallel for schedule(static) private(k, kp)
//#pragma omp parallel for
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
        for (k = 0; k < LDB; k++) {
            dataType Btmp[cnt];
            for (i = 0; i < cnt; i++) {
                Btmp[i] = B[swaped[index[i]] + k * LDB];
            }
            for (i = 0; i < cnt; i++) {
                B[index[i] + k * LDB] = Btmp[i];
            }
        }
        /*
         *        Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
         */
        TRSM_("L", "L", "N", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
        /*
         *        Compute D \ B -> B   [ D \ (L \P**T * B) ]
         */
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
            //#pragma omp parallel for schedule(static) private(bkm1, bk)
            for (j = 1; j <= NRHS; j++) {
                bkm1 = B[i - 1 + (j - 1) * LDB] / akm1k;
                bk = B[i + (j - 1) * LDB] / akm1k;
                B[i - 1 + (j - 1) * LDB] = (ak * bkm1 - bk) / deNom;
                B[i + (j - 1) * LDB] = (akm1 * bk - bkm1) / deNom;
            }
        }
        /*
         *        Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
         */
        TRSM_("L", "L", "T", "U", &N, &NRHS, &CONE, A, &LDA, B, &LDB);
/*
 *        P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
 *
 *        Interchange rows K and IPIV(K) of matrix B in reverse order
 *        from the formation order of IPIVi vector for Lower case.
 *
 *        (We can do the simple loop over IPIV with decrement -1,
 *        since the ABS value of IPIVi represents the row index
 *        of the interchange with row i in both 1x1 and 2x2 pivot cases)
 */
//#pragma omp parallel for schedule(static) private(k, kp)
//#pragma omp parallel for
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
        for (k = 0; k < LDB; k++) {
            dataType Btmp[cnt];
            for (i = 0; i < cnt; i++) {
                Btmp[i] = B[swaped[index[i]] + k * LDB];
            }
            for (i = 0; i < cnt; i++) {
                B[index[i] + k * LDB] = Btmp[i];
            }
        }
        /*
         *        END Lower
         */
    }
    return;
    /*
     *     End of SSYTRS_3
     */
}