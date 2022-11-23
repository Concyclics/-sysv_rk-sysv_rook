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

    int step = 32 / sizeof(dataType);
#if defined COMPLEX16
    if (N < 18000) {
        step *= 2;
    }
#endif
    int numThread = NRHS / step;
    if (numThread == 0) {
        step = 0;
        numThread = 1;
    }
    int last = (numThread - 1) * step;
    int length;

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
        Xerbla("SSYTRS_ROOK", &neg_info, 11);
#endif
#ifdef DOUBLE
        Xerbla("DSYTRS_ROOK", &neg_info, 11);
#endif
#ifdef COMPLEX
        Xerbla("CSYTRS_ROOK", &neg_info, 11);
#endif
#ifdef COMPLEX16
        Xerbla("ZSYTRS_ROOK", &neg_info, 11);
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
#pragma omp parallel for private(i, kp, akm1k, akm1, ak, deNom, bkm1, bk, \
                                 intTmp, length, k)
        for (i = 0; i < numThread; i++) {
            k = N;
            if (i < numThread - 1) {
                length = step;
            } else {
                length = NRHS - last;
            }
            while (k > 0) {
                if (ipiv[k - 1] > 0) {
                    kp = ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + k - 1 + i * step * LDB, ldb,
                              B + kp - 1 + i * step * LDB, ldb);
                    }
                    intTmp = k - 1;
                    GER_(&intTmp, &length, &NEG_CONE, A + (k - 1) * LDA,
                         &intOne, B + k - 1 + i * step * LDB, ldb,
                         B + i * step * LDB, ldb);
                    dataType ALPHA_ = T_ONE / A[k - 1 + (k - 1) * (LDA)];
#if defined(COMPLEX) || defined(COMPLEX16)
                    SCAL_(&length, (void*)&ALPHA_, B + k - 1 + i * step * LDB,
                          ldb);
#else
                    SCAL_(&length, &ALPHA_, B + k - 1 + i * step * LDB, ldb);
#endif

                    k -= 1;
                } else {
                    kp = -ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    kp = -ipiv[k - 2];
                    if (kp != k - 1) {
                        SWAP_(&length, B + (k - 2) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    if (k > 2) {
                        intTmp = k - 2;
                        GER_(&intTmp, &length, &NEG_CONE, A + (k - 1) * LDA,
                             &intOne, B + (k - 1) + i * step * LDB, &LDB,
                             B + i * step * LDB, &LDB);
                        GER_(&intTmp, &length, &NEG_CONE, A + (k - 2) * LDA,
                             &intOne, B + (k - 2) + i * step * LDB, &LDB,
                             B + i * step * LDB, &LDB);
                    }
                    akm1k = A[k - 2 + (k - 1) * (LDA)];
                    akm1 = A[k - 2 + (k - 2) * (LDA)] / akm1k;
                    ak = A[k - 1 + (k - 1) * (LDA)] / akm1k;
                    deNom = akm1 * ak - CONE;
                    for (j = i * step; j < i * step + length; j++) {
                        bkm1 = B[k - 2 + j * (LDB)] / akm1k;
                        bk = B[k - 1 + j * (LDB)] / akm1k;
                        B[k - 2 + j * (LDB)] = (bkm1 * ak - bk) / deNom;
                        B[k - 1 + j * (LDB)] = (akm1 * bk - bkm1) / deNom;
                    }
                    k -= 2;
                }
            }
        }

#pragma omp parallel for private(i, kp, intTmp, length, k)
        for (i = 0; i < numThread; i++) {
            k = 1;
            if (i < numThread - 1) {
                length = step;
            } else {
                length = NRHS - last;
            }
            while (k <= N) {
                if (ipiv[k - 1] > 0) {
                    if (k > 1) {
                        intTmp = k - 1;
                        GEMV_("T", &intTmp, &length, &NEG_CONE,
                              B + i * step * LDB, &LDB, A + (k - 1) * LDA,
                              &intOne, &CONE, B + k - 1 + i * step * LDB, &LDB);
                    }

                    kp = ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    k += 1;
                } else {
                    if (k > 1) {
                        intTmp = k - 1;
                        GEMV_("T", &intTmp, &length, &NEG_CONE,
                              B + i * step * LDB, &LDB, A + (k - 1) * LDA,
                              &intOne, &CONE, B + k - 1 + i * step * LDB, &LDB);
                        GEMV_("T", &intTmp, &length, &NEG_CONE,
                              B + i * step * LDB, &LDB, A + k * LDA, &intOne,
                              &CONE, B + k + i * step * LDB, &LDB);
                    }

                    kp = -ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    kp = -ipiv[k];
                    if (kp != k + 1) {
                        SWAP_(&length, B + k + i * step * LDB, &LDB,
                              B + kp - 1 + i * step * LDB, &LDB);
                    }
                    k += 2;
                }
            }
        }

    } else {
#pragma omp parallel for private(i, kp, akm1k, akm1, ak, deNom, bkm1, bk, \
                                 intTmp, length, k)
        for (i = 0; i < numThread; i++) {
            if (i < numThread - 1) {
                length = step;
            } else {
                length = NRHS - last;
            }
            k = 1;
            while (k <= N) {
                if (ipiv[k - 1] > 0) {
                    kp = ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }

                    if (k < N) {
                        intTmp = N - k;
                        GER_(&intTmp, &length, &NEG_CONE, A + k + (k - 1) * LDA,
                             &intOne, B + k - 1 + i * step * LDB, &LDB,
                             B + k + i * step * LDB, &LDB);
                    }

                    dataType ALPHA_ = T_ONE / A[k - 1 + (k - 1) * (LDA)];
#if defined(COMPLEX) || defined(COMPLEX16)
                    SCAL_(&length, (void*)&ALPHA_, B + k - 1 + i * step * LDB,
                          &LDB);
#else
                    SCAL_(&length, &ALPHA_, B + k - 1 + i * step * LDB, &LDB);
#endif
                    k += 1;
                } else {
                    kp = -ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    kp = -ipiv[k];
                    if (kp != k + 1) {
                        SWAP_(&length, B + k + i * step * LDB, &LDB,
                              B + kp - 1 + i * step * LDB, &LDB);
                    }

                    if (k < N - 1) {
                        intTmp = N - k - 1;
                        GER_(&intTmp, &length, &NEG_CONE,
                             A + k + 1 + (k - 1) * LDA, &intOne,
                             B + k - 1 + i * step * LDB, &LDB,
                             B + k + 1 + i * step * LDB, &LDB);
                        GER_(&intTmp, &length, &NEG_CONE, A + k + 1 + k * LDA,
                             &intOne, B + k + i * step * LDB, &LDB,
                             B + k + 1 + i * step * LDB, &LDB);
                    }

                    akm1k = A[k + (k - 1) * (LDA)];
                    akm1 = A[k - 1 + (k - 1) * (LDA)] / akm1k;
                    ak = A[k + k * (LDA)] / akm1k;
                    deNom = akm1 * ak - T_ONE;
                    for (j = i * step; j < i * step + length; j++) {
                        bkm1 = B[k - 1 + j * (LDB)] / akm1k;
                        bk = B[k + j * (LDB)] / akm1k;
                        B[k - 1 + j * (LDB)] = (bkm1 * ak - bk) / deNom;
                        B[k + j * (LDB)] = (bk * akm1 - bkm1) / deNom;
                    }
                    k += 2;
                }
            }
        }

#pragma omp parallel for private(i, kp, intTmp, length, k)
        for (i = 0; i < numThread; i++) {
            if (i < numThread - 1) {
                length = step;
            } else {
                length = NRHS - last;
            }
            k = N;
            while (k > 0) {
                if (ipiv[k - 1] > 0) {
                    if (k < N) {
                        intTmp = N - k;
                        GEMV_("T", &intTmp, &length, &NEG_CONE,
                              B + k + i * step * LDB, &LDB,
                              A + k + (k - 1) * LDA, &intOne, &CONE,
                              B + k - 1 + i * step * LDB, &LDB);
                    }

                    kp = ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    k -= 1;
                } else {
                    if (k < N) {
                        intTmp = N - k;
                        GEMV_("T", &intTmp, &length, &NEG_CONE,
                              B + k + i * step * LDB, &LDB,
                              A + k + (k - 1) * LDA, &intOne, &CONE,
                              B + k - 1 + i * step * LDB, &LDB);
                        GEMV_("T", &intTmp, &length, &NEG_CONE,
                              B + k + i * step * LDB, &LDB,
                              A + k + (k - 2) * LDA, &intOne, &CONE,
                              B + k - 2 + i * step * LDB, &LDB);
                    }
                    kp = -ipiv[k - 1];
                    if (kp != k) {
                        SWAP_(&length, B + (k - 1) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    kp = -ipiv[k - 2];
                    if (kp != k - 1) {
                        SWAP_(&length, B + (k - 2) + i * step * LDB, &LDB,
                              B + (kp - 1) + i * step * LDB, &LDB);
                    }
                    k -= 2;
                }
            }
        }
    }
}
