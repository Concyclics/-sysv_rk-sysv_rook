/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of SYTF2_RK.
 * Author: CHEN Han
 * Create: 2022-06-26
 *******************************************************************************/

#include "SYTF2_RK.h"

void SYTF2_RK(const char* uplo,
              const int* n,
              dataType* A,
              const int* lda,
              dataType* E,
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
    int upper, done;
    int i, imax, j, jmax, itemp, k, kk, kp, kstep, p, ii;
    dataAccu absakk, colmax, rowmax, stemp, sfmin;
    dataType d11, d12, d21, d22, t, wk, wkm1, wkp1, z;

    const int N = *n;
    const int LDA = *lda;

    int lsamax;
    dataAccu slamch;

    /* test the input parameters */

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
        int negInfo = -*info;
#ifdef SINGLE
        Xerbla("SSYTF2_RK", &negInfo, 7);
#endif
#ifdef DOUBLE
        Xerbla("DSYTF2_RK", &negInfo, 7);
#endif
#ifdef COMPLEX
        Xerbla("CSYTF2_RK", &negInfo, 7);
#endif
#ifdef COMPLEX16
        Xerbla("ZSYTF2_RK", &negInfo, 7);
#endif
        return;
    }

    dataAccu ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;

    sfmin = lamch_("S");
    if (upper) {
        E[0] = ZERO;

        k = N;

        while (k > 0) {
            kstep = 1;
            p = k;

            absakk = ABS_(A[k - 1 + (k - 1) * LDA]);

            if (k > 1) {
                intTmp = k - 1;
                imax = I_AMAX(&intTmp, A + (k - 1) * LDA, &intOne);
                colmax = ABS_(A[imax - 1 + (k - 1) * LDA]);
            } else {
                colmax = ZERO;
            }

            if (MAX(absakk, colmax) == ZERO) {

                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                if (k > 1) {
                    E[k - 1] = ZERO;
                }
            } else {

                if (absakk >= ALPHA * colmax) {
                    kp = k;
                } else {
                    done = false;

                    while (!done) {
                        if (imax != k) {
                            intTmp = k - imax;
                            jmax =
                                imax + I_AMAX(&intTmp,
                                              A + imax - 1 + (imax)*LDA, &LDA);
                            rowmax = ABS_(A[imax - 1 + (jmax - 1) * LDA]);
                        } else {
                            rowmax = ZERO;
                        }
                        if (imax > 1) {
                            intTmp = imax - 1;
                            itemp =
                                I_AMAX(&intTmp, A + (imax - 1) * LDA, &intOne);
                            stemp = ABS_(A[itemp - 1 + (imax - 1) * LDA]);
                            if (stemp > rowmax) {
                                rowmax = stemp;
                                jmax = itemp;
                            }
                        }

                        if ((ABS_(A[imax - 1 + (imax - 1) * LDA]) >=
                             ALPHA * rowmax)) {

                            kp = imax;
                            done = true;
                        }

                        else if (p == jmax || rowmax <= colmax) {

                            kp = imax;
                            done = true;
                            kstep = 2;
                        } else {
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                        }
                    }
                }
                if (kstep == 2 && p != k) {

                    if (p > 1) {
                        intTmp = p - 1;
                        SWAP_(&intTmp, A + (k - 1) * LDA, &intOne,
                              A + (p - 1) * LDA, &intOne);
                    }
                    if (p < k - 1) {
                        intTmp = k - p - 1;
                        SWAP_(&intTmp, A + p + (k - 1) * LDA, &intOne,
                              A + p - 1 + p * LDA, &LDA);
                    }
                    t = A[k - 1 + (k - 1) * LDA];
                    A[k - 1 + (k - 1) * LDA] = A[p - 1 + (p - 1) * LDA];
                    A[p - 1 + (p - 1) * LDA] = t;
                    if (k < N) {
                        intTmp = N - k;
                        SWAP_(&intTmp, A + k - 1 + k * LDA, &LDA,
                              A + p - 1 + k * LDA, &LDA);
                    }
                }

                kk = k - kstep + 1;
                if (kp != kk) {

                    if (kp > 1) {
                        intTmp = kp - 1;
                        SWAP_(&intTmp, A + (kk - 1) * LDA, &intOne,
                              A + (kp - 1) * LDA, &intOne);
                    }
                    if (kk > 1 && kp < kk - 1) {
                        intTmp = kk - kp - 1;
                        SWAP_(&intTmp, A + kp + (kk - 1) * LDA, &intOne,
                              A + kp - 1 + kp * LDA, &LDA);
                    }
                    t = A[kk - 1 + (kk - 1) * LDA];
                    A[kk - 1 + (kk - 1) * LDA] = A[kp - 1 + (kp - 1) * LDA];
                    A[kp - 1 + (kp - 1) * LDA] = t;
                    if (kstep == 2) {
                        t = A[k - 2 + (k - 1) * LDA];
                        A[k - 2 + (k - 1) * LDA] = A[kp - 1 + (k - 1) * LDA];
                        A[kp - 1 + (k - 1) * LDA] = t;
                    }
                    if (k < N) {
                        intTmp = N - k;
                        SWAP_(&intTmp, A + kk - 1 + k * LDA, &LDA,
                              A + kp - 1 + k * LDA, &LDA);
                    }
                }
                if (kstep == 1) {

                    if (k > 1) {
                        if (ABS_(A[k - 1 + (k - 1) * LDA]) >= sfmin) {
                            d11 = T_ONE / A[k - 1 + (k - 1) * LDA];
                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * LDA, &intOne,
                                 A, &LDA);
                            d11 = -d11;
                            intTmp = k - 1;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + (k - 1) * LDA, &intOne);
#endif
                        } else {
                            d11 = A[k - 1 + (k - 1) * LDA];

#pragma omp parallel for private(ii)
                            for (ii = 1; ii < k; ii++) {
                                A[ii - 1 + (k - 1) * LDA] /= d11;
                            }

                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * LDA, &intOne,
                                 A, &LDA);
                            d11 = -d11;
                        }

                        E[k - 1] = ZERO;
                    }
                } else {

                    if (k > 2) {
                        d12 = A[k - 2 + (k - 1) * LDA];
                        d22 = A[k - 2 + (k - 2) * LDA] / d12;
                        d11 = A[k - 1 + (k - 1) * LDA] / d12;
                        t = T_ONE / (d11 * d22 - T_ONE);
                        for (j = k - 2; j >= 1; --j) {
                            wkm1 = t * (d11 * A[j - 1 + (k - 2) * LDA] -
                                        A[j - 1 + (k - 1) * LDA]);
                            wk = t * (d22 * A[j - 1 + (k - 1) * LDA] -
                                      A[j - 1 + (k - 2) * LDA]);
#pragma omp parallel for private(i)
                            for (i = j; i >= 1; --i) {
                                A[i - 1 + (j - 1) * LDA] -=
                                    A[i - 1 + (k - 1) * LDA] / d12 * wk +
                                    A[i - 1 + (k - 2) * LDA] / d12 * wkm1;
                            }

                            A[j - 1 + (k - 1) * LDA] = wk / d12;
                            A[j - 1 + (k - 2) * LDA] = wkm1 / d12;
                        }
                    }
                    E[k - 1] = A[k - 2 + (k - 1) * LDA];
                    E[k - 2] = ZERO;
                    A[k - 2 + (k - 1) * LDA] = ZERO;
                }
            }
            if (kstep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k - 2] = -kp;
            }
            k = k - kstep;
        }
    } else {

        E[N - 1] = ZERO;

        k = 1;

        while (k <= N) {
            kstep = 1;
            p = k;

            absakk = ABS_(A[k - 1 + (k - 1) * LDA]);

            if (k < N) {
                intTmp = N - k;
                imax = k + I_AMAX(&intTmp, A + k + (k - 1) * LDA, &intOne);
                colmax = ABS_(A[imax - 1 + (k - 1) * LDA]);
            } else {
                colmax = ZERO;
            }

            if (MAX(absakk, colmax) == ZERO) {

                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                if (k < N) {
                    E[k - 1] = ZERO;
                }
            } else {

                if (absakk >= ALPHA * colmax) {

                    kp = k;
                } else {
                    done = false;
                    while (!done) {
                        if (imax != k) {
                            intTmp = imax - k;
                            jmax = k - 1 +
                                   I_AMAX(&intTmp, A + imax - 1 + (k - 1) * LDA,
                                          &LDA);
                            rowmax = ABS_(A[imax - 1 + (jmax - 1) * LDA]);
                        } else {
                            rowmax = ZERO;
                        }

                        if (imax < N) {
                            intTmp = N - imax;
                            itemp = imax + I_AMAX(&intTmp,
                                                  A + imax + (imax - 1) * LDA,
                                                  &intOne);
                            stemp = ABS_(A[itemp - 1 + (imax - 1) * LDA]);
                            if (stemp > rowmax) {
                                rowmax = stemp;
                                jmax = itemp;
                            }
                        }

                        if (ABS_(A[imax - 1 + (imax - 1) * LDA]) >=
                            ALPHA * rowmax) {

                            kp = imax;
                            done = true;
                        } else if (p == jmax || rowmax <= colmax) {

                            kp = imax;
                            kstep = 2;
                            done = true;
                        } else {
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                        }
                    }
                }

                if (kstep == 2 && p != k) {

                    if (p < N) {
                        intTmp = N - p;
                        SWAP_(&intTmp, A + p + (k - 1) * LDA, &intOne,
                              A + p + (p - 1) * LDA, &intOne);
                    }
                    if (p > k + 1) {
                        intTmp = p - k - 1;
                        SWAP_(&intTmp, A + k + (k - 1) * LDA, &intOne,
                              A + p - 1 + k * LDA, &LDA);
                    }
                    t = A[k - 1 + (k - 1) * LDA];
                    A[k - 1 + (k - 1) * LDA] = A[p - 1 + (p - 1) * LDA];
                    A[p - 1 + (p - 1) * LDA] = t;

                    if (k > 1) {
                        intTmp = k - 1;
                        SWAP_(&intTmp, A + k - 1, &LDA, A + p - 1, &LDA);
                    }
                }

                kk = k + kstep - 1;
                if (kp != kk) {

                    if (kp < N) {
                        intTmp = N - kp;
                        SWAP_(&intTmp, A + kp + (kk - 1) * LDA, &intOne,
                              A + kp + (kp - 1) * LDA, &intOne);
                    }
                    if (kk < N && kp > kk + 1) {
                        intTmp = kp - kk - 1;
                        SWAP_(&intTmp, A + kk + (kk - 1) * LDA, &intOne,
                              A + kp - 1 + kk * LDA, &LDA);
                    }
                    t = A[kk - 1 + (kk - 1) * LDA];
                    A[kk - 1 + (kk - 1) * LDA] = A[kp - 1 + (kp - 1) * LDA];
                    A[kp - 1 + (kp - 1) * LDA] = t;
                    if (kstep == 2) {
                        t = A[k + (k - 1) * LDA];
                        A[k + (k - 1) * LDA] = A[kp - 1 + (k - 1) * LDA];
                        A[kp - 1 + (k - 1) * LDA] = t;
                    }

                    if (k > 1) {
                        intTmp = k - 1;
                        SWAP_(&intTmp, A + kk - 1, &LDA, A + kp - 1, &LDA);
                    }
                }

                if (kstep == 1) {
                    if (k < N) {

                        if (ABS_(A[k - 1 + (k - 1) * LDA]) >= sfmin) {

                            d11 = T_ONE / A[k - 1 + (k - 1) * LDA];
                            intTmp = N - k;
                            d11 = -d11;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * LDA,
                                 &intOne, A + k + k * LDA, &LDA);
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
                            d11 = A[k - 1 + (k - 1) * LDA];

#pragma omp parallel for private(ii)
                            for (ii = k + 1; ii <= N; ii++) {
                                A[ii - 1 + (k - 1) * LDA] /= d11;
                            }
                            d11 = -d11;
                            intTmp = N - k;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * LDA,
                                 &intOne, A + k + k * LDA, &LDA);
                            d11 = -d11;
                        }
                        E[k - 1] = ZERO;
                    }
                } else {

                    if (k < N - 1) {
                        d21 = A[k + (k - 1) * LDA];
                        d11 = A[k + k * LDA] / d21;
                        d22 = A[k - 1 + (k - 1) * LDA] / d21;
                        t = T_ONE / (d11 * d22 - T_ONE);

                        for (j = k + 2; j <= N; ++j) {
                            wk = t * (d11 * A[j - 1 + (k - 1) * LDA] -
                                      A[j - 1 + k * LDA]);
                            wkp1 = t * (d22 * A[j - 1 + k * LDA] -
                                        A[j - 1 + (k - 1) * LDA]);

#pragma omp parallel for private(i)
                            for (i = j; i <= N; ++i) {
                                A[i - 1 + (j - 1) * LDA] -=
                                    wk * A[i - 1 + (k - 1) * LDA] / d21 +
                                    wkp1 * A[i - 1 + k * LDA] / d21;
                            }

                            A[j - 1 + (k - 1) * LDA] = wk / d21;
                            A[j - 1 + k * LDA] = wkp1 / d21;
                        }
                    }

                    E[k - 1] = A[k + (k - 1) * LDA];
                    E[k] = ZERO;
                    A[k + (k - 1) * LDA] = ZERO;
                }
            }
            if (kstep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k] = -kp;
            }
            k += kstep;
        }
    }
    return;
}
