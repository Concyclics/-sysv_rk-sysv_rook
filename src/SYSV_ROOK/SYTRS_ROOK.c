/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYTRS_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-27
 *******************************************************************************/

#include "SYTRS_ROOK.h"
void SYTRS_ROOK(const char* uplo,
                const int* n,
                const int* nrhs,
                dataType* A,
                const int* lda,
                int* ipiv,
                dataType* B,
                const int* ldb,
                int* info) {
    const dataType CONE = 1;
    const dataType NEG_CONE = -1;
    const int intOne = 1;
    int intTmp;

    const int N = *n;
    const int NRHS = *nrhs;
    const int LDA = *lda;
    const int LDB = *ldb;
    int upper;
    int i, j, k, kp;
    dataType ak, akm1, akm1k, bk, bkm1, deNom;

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
        *info = -8;
    }
    int neg_info = -*info;
    if (*info != 0) {
#ifdef SINGLE
        Xerbla("SYTRS_ROOK", &neg_info, 10);
#endif
#ifdef DOUBLE
        Xerbla("DYTRS_ROOK", &neg_info, 10);
#endif
#ifdef COMPLEX
        Xerbla("CYTRS_ROOK", &neg_info, 10);
#endif
#ifdef COMPLEX16
        Xerbla("ZYTRS_ROOK", &neg_info, 10);
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
        /*
         *  Solve A*X = B, where A = U*D*U**T.
         *
         *  First solve U*D*X = B, overwriting B with X.
         *
         *  K is the main loop index, decreasing from N to 1 in steps of
         *  1 or 2, depending on the size of the diagonal blocks.
         */
        k = N;
        while (k > 0) {
            if (ipiv[k - 1] > 0) {
                /*
                 *  1 x 1 diagonal block
                 *
                 *  Interchange rows K and IPIV(K).
                 */
                kp = ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                /*
                 *  Multiply by inv(U(K)), where U(K) is the transformation
                 *  stored in column K of A.
                 */
                intTmp = k - 1;
                GER_(&intTmp, &NRHS, &NEG_CONE, A + (k - 1) * LDA, &intOne,
                     B + (k - 1), &LDB, B, &LDB);
                /*
                 *  Multiply by the inverse of the diagonal block.
                 */
                dataType ALPHA_ = T_ONE / A[k - 1 + (k - 1) * (LDA)];
#if defined(COMPLEX) || defined(COMPLEX16)
                SCAL_(&NRHS, (void*)&ALPHA_, B + k - 1, &LDB);
#else
                SCAL_(&NRHS, &ALPHA_, B + k - 1, &LDB);
#endif
                k -= 1;
            } else {
                /*
                 *  2 x 2 diagonal block
                 *
                 *  Interchange rows K and -IPIV(K) THEN K-1 and -IPIV(K-1)
                 */
                kp = -ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                kp = -ipiv[k - 2];
                if (kp != k - 1) {
                    SWAP_(&NRHS, B + (k - 2), &LDB, B + (kp - 1), &LDB);
                }
                /*
                 *  Multiply by inv(U(K)), where U(K) is the transformation
                 *  stored in columns K-1 and K of A.
                 */
                if (k > 2) {
                    intTmp = k - 2;
                    GER_(&intTmp, &NRHS, &NEG_CONE, A + (k - 1) * LDA, &intOne,
                         B + (k - 1), &LDB, B, &LDB);
                    GER_(&intTmp, &NRHS, &NEG_CONE, A + (k - 2) * LDA, &intOne,
                         B + (k - 2), &LDB, B, &LDB);
                }
                /*
                 *  Multiply by the inverse of the diagonal block.
                 */
                akm1k = A[k - 2 + (k - 1) * (LDA)];
                akm1 = A[k - 2 + (k - 2) * (LDA)] / akm1k;
                ak = A[k - 1 + (k - 1) * (LDA)] / akm1k;
                deNom = akm1 * ak - CONE;
                for (j = 1; j <= NRHS; j++) {
                    bkm1 = B[k - 2 + (j - 1) * (LDB)] / akm1k;
                    bk = B[k - 1 + (j - 1) * (LDB)] / akm1k;
                    B[k - 2 + (j - 1) * (LDB)] = (bkm1 * ak - bk) / deNom;
                    B[k - 1 + (j - 1) * (LDB)] = (akm1 * bk - bkm1) / deNom;
                }
                k -= 2;
            }
        }
        /*
         *  Next solve U**T *X = B, overwriting B with X.
         *
         *  K is the main loop index, increasing from 1 to N in steps of
         *  1 or 2, depending on the size of the diagonal blocks.
         */
        k = 1;

        while (k <= N) {
            if (ipiv[k - 1] > 0) {
                /*
                 *  1 x 1 diagonal block
                 *
                 *  Multiply by inv(U(K)), where U(K) is the transformation
                 *  stored in column K of A.
                 */
                if (k > 1) {
                    intTmp = k - 1;
                    GEMV_("T", &intTmp, &NRHS, &NEG_CONE, B, &LDB,
                          A + (k - 1) * LDA, &intOne, &CONE, B + k - 1, &LDB);
                }
                /*
                 *  Interchange rows K and IPIV(K).
                 */
                kp = ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                k += 1;
            } else {
                /*
                 *  2 x 2 diagonal block
                 *
                 *  Multiply by inv(U**T(K+1)), where U(K+1) is the
                 * transformation stored in columns K and K+1 of A.
                 */
                if (k > 1) {
                    intTmp = k - 1;
                    GEMV_("T", &intTmp, &NRHS, &NEG_CONE, B, &LDB,
                          A + (k - 1) * LDA, &intOne, &CONE, B + k - 1, &LDB);
                    GEMV_("T", &intTmp, &NRHS, &NEG_CONE, B, &LDB, A + k * LDA,
                          &intOne, &CONE, B + k, &LDB);
                }
                /*
                 *  Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1).
                 */
                kp = -ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                kp = -ipiv[k];
                if (kp != k + 1) {
                    SWAP_(&NRHS, B + k, &LDB, B + kp - 1, &LDB);
                }
                k += 2;
            }
        }
    } else {
        /*
         *  Solve U**T *X = B, overwriting B with X.
         *
         *  K is the main loop index, increasing from N to 1 in steps of
         *  -1 or -2, depending on the size of the diagonal blocks.
         */
        k = 1;
        while (k <= N) {
            if (ipiv[k - 1] > 0) {
                /*
                 *  1 x 1 diagonal block
                 *
                 *  Interchange rows K and IPIV(K).
                 */
                kp = ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                /*
                 *  Multiply by inv(L(K)), where L(K) is the transformation
                 *  stored in column K of A.
                 */
                if (k < N) {
                    intTmp = N - k;
                    GER_(&intTmp, &NRHS, &NEG_CONE, A + k + (k - 1) * LDA,
                         &intOne, B + k - 1, &LDB, B + k, &LDB);
                }
                /*
                 *  Multiply by the inverse of the diagonal block.
                 */
                dataType ALPHA_ = T_ONE / A[k - 1 + (k - 1) * (LDA)];
#if defined(COMPLEX) || defined(COMPLEX16)
                SCAL_(&NRHS, (void*)&ALPHA_, B + k - 1, &LDB);
#else
                SCAL_(&NRHS, &ALPHA_, B + k - 1, &LDB);
#endif
                k += 1;
            } else {
                /*
                 *  2 x 2 diagonal block
                 *
                 *  Interchange rows K and -IPIV(K) THEN K+1 and -IPIV(K+1)
                 */
                kp = -ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                kp = -ipiv[k];
                if (kp != k + 1) {
                    SWAP_(&NRHS, B + k, &LDB, B + kp - 1, &LDB);
                }
                /*
                 *  Multiply by inv(L(K+1)), where L(K+1) is the transformation
                 *  stored in columns K and K+1 of A.
                 */
                if (k < N - 1) {
                    intTmp = N - k - 1;
                    GER_(&intTmp, &NRHS, &NEG_CONE, A + k + 1 + (k - 1) * LDA,
                         &intOne, B + k - 1, &LDB, B + k + 1, &LDB);
                    GER_(&intTmp, &NRHS, &NEG_CONE, A + k + 1 + k * LDA,
                         &intOne, B + k, &LDB, B + k + 1, &LDB);
                }
                /*
                 *  Multiply by the inverse of the diagonal block.
                 */
                akm1k = A[k + (k - 1) * (LDA)];
                akm1 = A[k - 1 + (k - 1) * (LDA)] / akm1k;
                ak = A[k + k * (LDA)] / akm1k;
                deNom = akm1 * ak - T_ONE;
                for (j = 1; j <= NRHS; j++) {
                    bkm1 = B[k - 1 + (j - 1) * (LDB)] / akm1k;
                    bk = B[k + (j - 1) * (LDB)] / akm1k;
                    B[k - 1 + (j - 1) * (LDB)] = (bkm1 * ak - bk) / deNom;
                    B[k + (j - 1) * (LDB)] = (bk * akm1 - bkm1) / deNom;
                }
                k += 2;
            }
        }
        /*
         *  Next solve L**T *X = B, overwriting B with X.
         *
         *  K is the main loop index, decreasing from N to 1 in steps of
         *  1 or 2, depending on the size of the diagonal blocks.
         */
        k = N;

        while (k > 0) {
            if (ipiv[k - 1] > 0) {
                /*
                 *  1 x 1 diagonal block
                 *
                 *  Multiply by inv(L(K)), where L(K) is the transformation
                 *  stored in column K of A.
                 */
                if (k < N) {
                    intTmp = N - k;
                    GEMV_("T", &intTmp, &NRHS, &NEG_CONE, B + k, &LDB,
                          A + k + (k - 1) * LDA, &intOne, &CONE, B + k - 1,
                          &LDB);
                }
                /*
                 *  Interchange rows K and IPIV(K).
                 */
                kp = ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                k -= 1;
            } else {
                /*
                 *  2 x 2 diagonal block
                 *
                 *  Multiply by inv(L**T(K-1)), where L(K-1) is the
                 * transformation stored in columns K-1 and K of A.
                 */
                if (k < N) {
                    intTmp = N - k;
                    GEMV_("T", &intTmp, &NRHS, &NEG_CONE, B + k, &LDB,
                          A + k + (k - 1) * LDA, &intOne, &CONE, B + k - 1,
                          &LDB);
                    GEMV_("T", &intTmp, &NRHS, &NEG_CONE, B + k, &LDB,
                          A + k + (k - 2) * LDA, &intOne, &CONE, B + k - 2,
                          &LDB);
                }
                kp = -ipiv[k - 1];
                if (kp != k) {
                    SWAP_(&NRHS, B + (k - 1), &LDB, B + (kp - 1), &LDB);
                }
                kp = -ipiv[k - 2];
                if (kp != k - 1) {
                    SWAP_(&NRHS, B + (k - 2), &LDB, B + (kp - 1), &LDB);
                }
                k -= 2;
            }
        }
    }
    // End of SYTRS_ROOK
}