/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYTF2_ROOK.
 * Author: CHEN Han
 * Create: 2022-06-27
 *******************************************************************************/

#include "SYTF2_ROOK.h"

void SYTF2_ROOK(const char* uplo,
                const int* n,
                dataType* A,
                const int* lda,
                int* ipiv,
                int* info) {
    const int intOne = 1;
    int intTmp;
#ifdef COMPLEX
    const dataType CONE = 1;
    const dataType NEG_CONE = -1;
#endif
#ifdef COMPLEX16
    const dataType CONE = 1;
    const dataType NEG_CONE = -1;
#endif

    const int LDA = *lda;
    const int N = *n;

    int upper, done;
    int i, iMax, j, jMax, iTemp, k, kk, kp, kStep, p, ii;
    dataAccu absAkk, Alpha, colMax, rowMax, sTemp, sfMin;
    dataType d11, d12, d21, d22, t, wk, wkm1, wkp1, z;

    Alpha = (ONE + sqrt(SEVTEN)) / EIGHT;

    dataAccu slamch;

    *info = 0;

    upper = KmlLsame(*uplo, 'U');
    if (!upper && !KmlLsame(*uplo, 'L')) {
        *info = -1;
    } else if (N < 0) {
        *info = -2;
    } else if (LDA < MAX(1, N)) {
        *info = -4;
    }

    if (*info != 0) {
        int neg_info = -*info;
#ifdef SINGLE
        Xerbla("SSYTF2_ROOK", &neg_info, 11);
#endif
#ifdef DOUBLE
        Xerbla("DSYTF2_ROOK", &neg_info, 11);
#endif
#ifdef COMPLEX
        Xerbla("CSYTF2_ROOK", &neg_info, 11);
#endif
#ifdef COMPLEX16
        Xerbla("ZSYTF2_ROOK", &neg_info, 11);
#endif
        return;
    }

    sfMin = lamch_("S");
    if (upper) {
        k = N;

        while (k > 0) {
            kStep = 1;
            p = k;

            absAkk = ABS_(A[k - 1 + (k - 1) * (LDA)]);

            if (k > 1) {
                intTmp = k - 1;
                iMax = I_AMAX(&intTmp, A + (k - 1) * (LDA), &intOne);
                colMax = ABS_(A[iMax - 1 + (k - 1) * (LDA)]);
            } else {
                colMax = ZERO;
            }

            if (MAX(absAkk, colMax) == ZERO) {
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
            } else {
                if (absAkk >= Alpha * colMax) {
                    kp = k;
                } else {
                    done = false;
                    while (!done) {
                        if (iMax != k) {
                            intTmp = k - iMax;
                            jMax = iMax + I_AMAX(&intTmp,
                                                 A + iMax - 1 + (iMax) * (LDA),
                                                 &LDA);
                            rowMax = ABS_(A[iMax - 1 + (jMax - 1) * (LDA)]);
                        } else {
                            rowMax = ZERO;
                        }
                        if (iMax > 1) {
                            intTmp = iMax - 1;
                            iTemp = I_AMAX(&intTmp, A + (iMax - 1) * (LDA),
                                           &intOne);
                            sTemp = ABS_(A[iTemp - 1 + (iMax - 1) * (LDA)]);
                            if (sTemp > rowMax) {
                                rowMax = sTemp;
                                jMax = iTemp;
                            }
                        }

                        if ((ABS_(A[iMax - 1 + (iMax - 1) * (LDA)]) >=
                             Alpha * rowMax)) {
                            kp = iMax;
                            done = true;
                        }

                        else if (p == jMax || rowMax <= colMax) {
                            kp = iMax;
                            kStep = 2;
                            done = true;
                        } else {
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                        }
                    }
                }

                if (kStep == 2 && p != k) {
                    if (p > 1) {
                        intTmp = p - 1;
                        SWAP_(&intTmp, A + (k - 1) * (LDA), &intOne,
                              A + (p - 1) * (LDA), &intOne);
                    }
                    if (p < k - 1) {
                        intTmp = k - p - 1;
                        SWAP_(&intTmp, A + p + (k - 1) * (LDA), &intOne,
                              A + p - 1 + p * (LDA), &LDA);
                    }
                    t = A[k - 1 + (k - 1) * (LDA)];
                    A[k - 1 + (k - 1) * (LDA)] = A[p - 1 + (p - 1) * (LDA)];
                    A[p - 1 + (p - 1) * (LDA)] = t;
                }

                kk = k - kStep + 1;
                if (kp != kk) {
                    if (kp > 1) {
                        intTmp = kp - 1;
                        SWAP_(&intTmp, A + (kk - 1) * (LDA), &intOne,
                              A + (kp - 1) * LDA, &intOne);
                    }
                    if (kk > 1 && kp < kk - 1) {
                        intTmp = kk - kp - 1;
                        SWAP_(&intTmp, A + kp + (kk - 1) * (LDA), &intOne,
                              A + kp - 1 + kp * (LDA), &LDA);
                    }
                    t = A[kk - 1 + (kk - 1) * (LDA)];
                    A[kk - 1 + (kk - 1) * (LDA)] = A[kp - 1 + (kp - 1) * (LDA)];
                    A[kp - 1 + (kp - 1) * (LDA)] = t;
                    if (kStep == 2) {
                        t = A[k - 2 + (k - 1) * (LDA)];
                        A[k - 2 + (k - 1) * (LDA)] =
                            A[kp - 1 + (k - 1) * (LDA)];
                        A[kp - 1 + (k - 1) * (LDA)] = t;
                    }
                }

                if (kStep == 1) {
                    if (k > 1) {
                        if (ABS_(A[k - 1 + (k - 1) * (LDA)]) >= sfMin) {
                            d11 = T_ONE / A[k - 1 + (k - 1) * (LDA)];
                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * (LDA),
                                 &intOne, A, &LDA);
                            d11 = -d11;
                            intTmp = k - 1;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + (k - 1) * LDA, &intOne);
#endif
                        } else {
                            d11 = A[k - 1 + (k - 1) * (LDA)];
#pragma omp parallel for private(ii)
                            for (ii = 1; ii < k; ii++) {
                                A[ii - 1 + (k - 1) * (LDA)] /= d11;
                            }
                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * (LDA),
                                 &intOne, A, &LDA);
                            d11 = -d11;
                        }
                    }
                } else {
                    if (k > 2) {
                        d12 = A[k - 2 + (k - 1) * (LDA)];
                        d22 = A[k - 2 + (k - 2) * (LDA)] / d12;
                        d11 = A[k - 1 + (k - 1) * (LDA)] / d12;
                        t = T_ONE / (d11 * d22 - T_ONE);
                        for (j = k - 2; j >= 1; j--) {
                            wkm1 = t * (d11 * A[j - 1 + (k - 2) * (LDA)] -
                                        A[j - 1 + (k - 1) * (LDA)]);
                            wk = t * (d22 * A[j - 1 + (k - 1) * (LDA)] -
                                      A[j - 1 + (k - 2) * (LDA)]);
#pragma omp parallel for private(i)
                            for (i = j; i >= 1; i--) {
                                A[i - 1 + (j - 1) * (LDA)] -=
                                    A[i - 1 + (k - 1) * (LDA)] / d12 * wk +
                                    A[i - 1 + (k - 2) * (LDA)] / d12 * wkm1;
                            }
                            A[j - 1 + (k - 1) * (LDA)] = wk / d12;
                            A[j - 1 + (k - 2) * (LDA)] = wkm1 / d12;
                        }
                    }
                }
            }
            if (kStep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k - 2] = -kp;
            }
            k = k - kStep;
        }
    } else {
        k = 1;

        while (k <= N) {
            kStep = 1;
            p = k;
            absAkk = ABS_(A[k - 1 + (k - 1) * (LDA)]);
            if (k < N) {
                intTmp = N - k;
                iMax = k + I_AMAX(&intTmp, A + k + (k - 1) * (LDA), &intOne);
                colMax = ABS_(A[iMax - 1 + (k - 1) * (LDA)]);
            } else {
                colMax = ZERO;
            }

            if (MAX(absAkk, colMax) == ZERO) {
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
            } else {
                if (absAkk >= Alpha * colMax) {
                    kp = k;
                } else {
                    done = false;
                    while (!done) {
                        if (iMax != k) {
                            intTmp = iMax - k;
                            jMax = k - 1 +
                                   I_AMAX(&intTmp,
                                          A + iMax - 1 + (k - 1) * (LDA), &LDA);
                            rowMax = ABS_(A[iMax - 1 + (jMax - 1) * (LDA)]);
                        } else {
                            rowMax = ZERO;
                        }

                        if (iMax < N) {
                            intTmp = N - iMax;
                            iTemp = iMax + I_AMAX(&intTmp,
                                                  A + iMax + (iMax - 1) * (LDA),
                                                  &intOne);
                            sTemp = ABS_(A[iTemp - 1 + (iMax - 1) * (LDA)]);
                            if (sTemp > rowMax) {
                                rowMax = sTemp;
                                jMax = iTemp;
                            }
                        }
                        if (ABS_(A[iMax - 1 + (iMax - 1) * (LDA)]) >=
                            Alpha * rowMax) {
                            kp = iMax;
                            done = true;
                        } else if (p == jMax || rowMax <= colMax) {
                            kp = iMax;
                            kStep = 2;
                            done = true;
                        } else {
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                        }
                    }
                }
                if (kStep == 2 && p != k) {
                    if (p < N) {
                        intTmp = N - p;
                        SWAP_(&intTmp, A + p + (k - 1) * (LDA), &intOne,
                              A + p + (p - 1) * (LDA), &intOne);
                    }
                    if (p > k + 1) {
                        intTmp = p - k - 1;
                        SWAP_(&intTmp, A + k + (k - 1) * (LDA), &intOne,
                              A + p - 1 + (k) * (LDA), &LDA);
                    }
                    t = A[k - 1 + (k - 1) * (LDA)];
                    A[k - 1 + (k - 1) * (LDA)] = A[p - 1 + (p - 1) * (LDA)];
                    A[p - 1 + (p - 1) * (LDA)] = t;
                }

                kk = k + kStep - 1;
                if (kp != kk) {
                    if (kp < N) {
                        intTmp = N - kp;
                        SWAP_(&intTmp, A + kp + (kk - 1) * (LDA), &intOne,
                              A + kp + (kp - 1) * (LDA), &intOne);
                    }
                    if (kk < N && kp > kk + 1) {
                        intTmp = kp - kk - 1;
                        SWAP_(&intTmp, A + kk + (kk - 1) * (LDA), &intOne,
                              A + kp - 1 + (kk) * (LDA), &LDA);
                    }
                    t = A[kk - 1 + (kk - 1) * (LDA)];
                    A[kk - 1 + (kk - 1) * (LDA)] = A[kp - 1 + (kp - 1) * (LDA)];
                    A[kp - 1 + (kp - 1) * (LDA)] = t;
                    if (kStep == 2) {
                        t = A[k + (k - 1) * (LDA)];
                        A[k + (k - 1) * (LDA)] = A[kp - 1 + (k - 1) * (LDA)];
                        A[kp - 1 + (k - 1) * (LDA)] = t;
                    }
                }
                if (kStep == 1) {
                    if (k < N) {
                        if (ABS_(A[k - 1 + (k - 1) * (LDA)]) >= sfMin) {
                            d11 = T_ONE / A[k - 1 + (k - 1) * (LDA)];
                            d11 = -d11;
                            intTmp = N - k;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * (LDA),
                                 &intOne, A + k + (k) * (LDA), &LDA);
                            d11 = -d11;
                            intTmp = N - k;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + k + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + k + (k - 1) * LDA,
                                  &intOne);
#endif
                        } else {
                            d11 = A[k - 1 + (k - 1) * (LDA)];
#pragma omp parallel for private(ii)
                            for (ii = k + 1; ii <= N; ii++) {
                                A[ii - 1 + (k - 1) * (LDA)] /= d11;
                            }

                            d11 = -d11;
                            intTmp = N - k;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * (LDA),
                                 &intOne, A + k + (k) * (LDA), &LDA);
                            d11 = -d11;
                        }
                    }
                } else {
                    if (k < N - 1) {
                        d21 = A[k + (k - 1) * (LDA)];
                        d11 = A[k + (k) * (LDA)] / d21;
                        d22 = A[k - 1 + (k - 1) * (LDA)] / d21;
                        t = T_ONE / (d11 * d22 - T_ONE);
                        for (j = k + 2; j <= N; j++) {
                            wk = t * (d11 * A[j - 1 + (k - 1) * (LDA)] -
                                      A[j - 1 + (k) * (LDA)]);
                            wkp1 = t * (d22 * A[j - 1 + (k) * (LDA)] -
                                        A[j - 1 + (k - 1) * (LDA)]);

#pragma omp parallel for private(i)
                            for (i = j; i <= N; i++) {
                                A[i - 1 + (j - 1) * (LDA)] -=
                                    wk * A[i - 1 + (k - 1) * (LDA)] / d21 +
                                    wkp1 * A[i - 1 + (k) * (LDA)] / d21;
                            }
                            A[j - 1 + (k - 1) * (LDA)] = wk / d21;
                            A[j - 1 + (k) * (LDA)] = wkp1 / d21;
                        }
                    }
                }
            }
            if (kStep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k] = -kp;
            }
            k += kStep;
        }
    }
    return;
}
