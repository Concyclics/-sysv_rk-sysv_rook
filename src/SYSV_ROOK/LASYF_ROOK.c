/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of LASYF_RK.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#include "LASYF_ROOK.h"

void LASYF_ROOK(const char* uplo,
                const int* n,
                const int* nb,
                int* kb,
                dataType* a,
                const int* lda,
                int* ipiv,
                dataType* w,
                const int* ldw,
                int* info) {
    const int intOne = 1;
    int intTmp, Tmp2;
    const dataType CONE = 1;
    const dataType NEG_CONE = -1;
    const int N = *n;
    const int NB = *nb;
    const int LDA = *lda;
    const int LDW = *ldw;
    int k, kw, kStep, p, iMax, jMax, iTemp, kk, kp, kkw, jb, jp1, jp2;
    int j, ii, jj, flag;
    dataAccu sfMin, absAKK, colMax, rowMax, sTemp;
    dataType d11, d12, d21, d22, r1, t;
    *info = 0;
    const dataAccu ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;
    char s = 'S';
    sfMin = lamch_(&s);
    if (KmlLsame(*uplo, 'U')) {
        k = N;
        while (1) {
            kw = NB + k - N;
            if ((k <= N - NB + 1 && NB < N) || k < 1) {
                break;
            }
            kStep = 1;
            p = k;
            COPY_(&k, a + (k - 1) * LDA, &intOne, w + (kw - 1) * LDW, &intOne);
            if (k < N) {
                intTmp = N - k;
                GEMV_("N", &k, &intTmp, &NEG_CONE, a + k * LDA, &LDA,
                      w + k - 1 + kw * LDW, &LDW, &CONE, w + (kw - 1) * LDW,
                      &intOne);
            }

            absAKK = ABS_(w[k - 1 + (kw - 1) * LDW]);

            if (k > 1) {
                intTmp = k - 1;
                iMax = I_AMAX(&intTmp, w + (kw - 1) * LDW, &intOne);
                colMax = ABS_(w[iMax - 1 + (kw - 1) * LDW]);
            } else {
                colMax = ZERO;
            }
            if (MAX(absAKK, colMax) == ZERO) {
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                COPY_(&k, w + (kw - 1) * LDW, &intOne, a + (k - 1) * LDA,
                      &intOne);
            } else {

                if (absAKK >= ALPHA * colMax) {
                    kp = k;
                } else {
                    flag = 1;
                    while (flag) {
                        COPY_(&iMax, a + (iMax - 1) * LDA, &intOne,
                              w + (kw - 2) * LDW, &intOne);
                        intTmp = k - iMax;
                        COPY_(&intTmp, a + iMax - 1 + iMax * LDA, &LDA,
                              w + iMax + (kw - 2) * LDW, &intOne);
                        if (k < N) {
                            intTmp = N - k;
                            GEMV_("N", &k, &intTmp, &NEG_CONE, a + k * LDA,
                                  &LDA, w + iMax - 1 + kw * LDW, &LDW, &CONE,
                                  w + (kw - 2) * LDW, &intOne);
                        }

                        if (iMax != k) {
                            intTmp = k - iMax;
                            jMax = iMax + I_AMAX(&intTmp,
                                                 w + iMax + (kw - 2) * LDW,
                                                 &intOne);
                            rowMax = ABS_(w[jMax - 1 + (kw - 2) * LDW]);
                        } else {
                            rowMax = ZERO;
                        }

                        if (iMax > 1) {
                            intTmp = iMax - 1;
                            iTemp =
                                I_AMAX(&intTmp, w + (kw - 2) * LDW, &intOne);
                            sTemp = ABS_(w[iTemp - 1 + (kw - 2) * LDW]);
                            if (sTemp > rowMax) {
                                rowMax = sTemp;
                                jMax = iTemp;
                            }
                        }
                        if (ABS_(w[iMax - 1 + (kw - 2) * LDW]) >=
                            ALPHA * rowMax) {
                            kp = iMax;
                            COPY_(&k, w + (kw - 2) * LDW, &intOne,
                                  w + (kw - 1) * LDW, &intOne);
                            flag = 0;
                        } else if ((p == jMax) ||
                                   (rowMax <=
                                    colMax))
                        {
                            kp = iMax;
                            kStep = 2;
                            flag = 0;
                        } else {
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                            COPY_(&k, w + (kw - 2) * LDW, &intOne,
                                  w + (kw - 1) * LDW, &intOne);
                        }
                    }
                }

                kk = k - kStep + 1;
                kkw = NB + kk - N;
                if ((kStep == 2) && (p != k)) {
                    intTmp = k - p;
                    COPY_(&intTmp, a + p + (k - 1) * LDA, &intOne,
                          a + p - 1 + p * LDA, &LDA);
                    COPY_(&p, a + (k - 1) * LDA, &intOne, a + (p - 1) * LDA,
                          &intOne);
                    intTmp = N - k + 1;
                    SWAP_(&intTmp, a + k - 1 + (k - 1) * LDA, &LDA,
                          a + p - 1 + (k - 1) * LDA, &LDA);
                    intTmp = N - kk + 1;
                    SWAP_(&intTmp, w + k - 1 + (kkw - 1) * LDW, &LDW,
                          w + p - 1 + (kkw - 1) * LDW, &LDW);
                }

                if (kp != kk) {
                    a[kp - 1 + (k - 1) * LDA] = a[kk - 1 + (k - 1) * LDA];
                    intTmp = k - 1 - kp;
                    COPY_(&intTmp, a + kp + (kk - 1) * LDA, &intOne,
                          a + kp - 1 + kp * LDA, &LDA);
                    COPY_(&kp, a + (kk - 1) * LDA, &intOne, a + (kp - 1) * LDA,
                          &intOne);
                    intTmp = N - kk + 1;
                    SWAP_(&intTmp, a + kk - 1 + (kk - 1) * LDA, &LDA,
                          a + kp - 1 + (kk - 1) * LDA, &LDA);
                    intTmp = N - kk + 1;
                    SWAP_(&intTmp, w + kk - 1 + (kkw - 1) * LDW, &LDW,
                          w + kp - 1 + (kkw - 1) * LDW, &LDW);
                }

                if (kStep == 1) {
                    COPY_(&k, w + (kw - 1) * LDW, &intOne, a + (k - 1) * LDA,
                          &intOne);
                    if (k > 1) {
                        if (ABS_(a[k - 1 + (k - 1) * LDA]) >= sfMin) {
                            r1 = T_ONE / a[k - 1 + (k - 1) * LDA];
                            intTmp = k - 1;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&r1, a + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &r1, a + (k - 1) * LDA, &intOne);
#endif
                        } else if (a[k - 1 + (k - 1) * LDA] != ZERO) {
#pragma omp parallel for private(ii)
                            for (ii = 1; ii <= k - 1; ii++) {
                                a[ii - 1 + (k - 1) * LDA] =
                                    a[ii - 1 + (k - 1) * LDA] /
                                    a[k - 1 + (k - 1) * LDA];
                            }
                        }
                    }
                } else {

                    if (k > 2) {
                        d12 = w[k - 2 + (kw - 1) * LDW];
                        d11 = w[k - 1 + (kw - 1) * LDW] / d12;
                        d22 = w[k - 2 + (kw - 2) * LDW] / d12;
                        t = T_ONE / (d11 * d22 - T_ONE);
#pragma omp parallel for private(j)
                        for (int j = 1; j <= k - 2; j++) {
                            a[j - 1 + (k - 2) * LDA] =
                                t * ((d11 * w[j - 1 + (kw - 2) * LDW] -
                                      w[j - 1 + (kw - 1) * LDW]) /
                                     d12);
                            a[j - 1 + (k - 1) * LDA] =
                                t * ((d22 * w[j - 1 + (kw - 1) * LDW] -
                                      w[j - 1 + (kw - 2) * LDW]) /
                                     d12);
                        }
                    }

                    a[k - 2 + (k - 2) * LDA] = w[k - 2 + (kw - 2) * LDW];
                    a[k - 2 + (k - 1) * LDA] = w[k - 2 + (kw - 1) * LDW];
                    a[k - 1 + (k - 1) * LDA] = w[k - 1 + (kw - 1) * LDW];
                }
            }
            if (kStep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k - 2] = -kp;
            }
            k -= kStep;
        }

#pragma omp parallel for private(j, jj, jb, intTmp, Tmp2)
        for (j = ((k - 1) / NB) * NB + 1; j >= 1; j -= NB) {
            jb = MIN(NB, k - j + 1);
            for (jj = j; jj <= j + jb - 1; jj++) {
                intTmp = jj - j + 1;
                Tmp2 = N - k;
                GEMV_("N", &intTmp, &Tmp2, &NEG_CONE, a + j - 1 + k * LDA, &LDA,
                      w + jj - 1 + kw * LDW, &LDW, &CONE,
                      a + j - 1 + (jj - 1) * LDA, &intOne);
            }
            if (j >= 2) {
                intTmp = j - 1;
                Tmp2 = N - k;
                GEMM_("N", "T", &intTmp, &jb, &Tmp2, &NEG_CONE, a + k * LDA,
                      &LDA, w + j - 1 + kw * LDW, &LDW, &CONE,
                      a + (j - 1) * LDA, &LDA);
            }
        }

        j = k + 1;
        while (j <= N) {
            kStep = 1;
            jp1 = 1;
            jj = j;
            jp2 = ipiv[j - 1];
            if (jp2 < 0) {
                jp2 = -jp2;
                j++;
                jp1 = -ipiv[j - 1];
                kStep = 2;
            }
            j++;
            if (jp2 != jj && j <= N) {
                intTmp = N - j + 1;
                SWAP_(&intTmp, a + jp2 - 1 + (j - 1) * LDA, &LDA,
                      a + jj - 1 + (j - 1) * LDA, &LDA);
            }
            jj = j - 1;
            if (jp1 != jj && kStep == 2) {
                intTmp = N - j + 1;
                SWAP_(&intTmp, a + jp1 - 1 + (j - 1) * LDA, &LDA,
                      a + jj - 1 + (j - 1) * LDA, &LDA);
            }
        }
        *kb = N - k;
    } else {
        k = 1;
        while (1) {
            if ((k >= NB && NB < N) || k > N) {
                break;
            }
            kStep = 1;
            p = k;
            intTmp = N - k + 1;
            COPY_(&intTmp, a + k - 1 + (k - 1) * LDA, &intOne,
                  w + k - 1 + (k - 1) * LDW, &intOne);
            if (k > 1) {
                intTmp = N - k + 1;
                Tmp2 = k - 1;
                GEMV_("N", &intTmp, &Tmp2, &NEG_CONE, a + k - 1, &LDA,
                      w + k - 1, &LDW, &CONE, w + k - 1 + (k - 1) * LDW,
                      &intOne);
            }
            absAKK = ABS_(w[k - 1 + (k - 1) * LDW]);

            if (k < N) {
                intTmp = N - k;
                iMax = k + I_AMAX(&intTmp, w + k + (k - 1) * LDW, &intOne);
                colMax = ABS_(w[iMax - 1 + (k - 1) * LDW]);
            } else {
                colMax = ZERO;
            }
            if (MAX(absAKK, colMax) == ZERO) {
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                intTmp = N - k + 1;
                COPY_(&intTmp, w + k - 1 + (k - 1) * LDW, &intOne,
                      a + k - 1 + (k - 1) * LDA, &intOne);
            } else {
                if (absAKK >= ALPHA * colMax) {
                    kp = k;
                } else {
                    flag = 1;
                    while (flag) {
                        intTmp = iMax - k;
                        COPY_(&intTmp, a + iMax - 1 + (k - 1) * LDA, &LDA,
                              w + k - 1 + k * LDW, &intOne);
                        intTmp = N - iMax + 1;
                        COPY_(&intTmp, a + iMax - 1 + (iMax - 1) * LDA, &intOne,
                              w + iMax - 1 + k * LDW, &intOne);
                        if (k > 1) {
                            intTmp = N - k + 1;
                            Tmp2 = k - 1;
                            GEMV_("N", &intTmp, &Tmp2, &NEG_CONE, a + k - 1,
                                  &LDA, w + iMax - 1, &LDW, &CONE,
                                  w + k - 1 + k * LDW, &intOne);
                        }
                        if (iMax != k) {
                            intTmp = iMax - k;
                            jMax =
                                k - 1 +
                                I_AMAX(&intTmp, w + k - 1 + k * LDW, &intOne);
                            rowMax = ABS_(w[jMax - 1 + k * LDW]);
                        } else {
                            rowMax = ZERO;
                        }

                        if (iMax < N) {
                            intTmp = N - iMax;
                            iTemp = iMax + I_AMAX(&intTmp, w + iMax + k * LDW,
                                                  &intOne);
                            sTemp = ABS_(w[iTemp - 1 + k * LDW]);
                            if (sTemp > rowMax) {
                                rowMax = sTemp;
                                jMax = iTemp;
                            }
                        }
                        if (ABS_(w[iMax - 1 + k * LDW]) >= ALPHA * rowMax) {

                            kp = iMax;

                            intTmp = N - k + 1;
                            COPY_(&intTmp, w + k - 1 + k * LDW, &intOne,
                                  w + k - 1 + (k - 1) * LDW, &intOne);
                            flag = 0;
                        }

                        else if ((p == jMax) || (rowMax <= colMax)) {
                            kp = iMax;
                            kStep = 2;
                            flag = 0;
                        } else {
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                            intTmp = N - k + 1;
                            COPY_(&intTmp, w + k - 1 + k * LDW, &intOne,
                                  w + k - 1 + (k - 1) * LDW, &intOne);
                        }
                    }
                }
                kk = k + kStep - 1;
                if ((kStep == 2) && (p != k)) {
                    intTmp = p - k;
                    COPY_(&intTmp, a + k - 1 + (k - 1) * LDA, &intOne,
                          a + p - 1 + (k - 1) * LDA, &LDA);
                    intTmp = N - p + 1;
                    COPY_(&intTmp, a + p - 1 + (k - 1) * LDA, &intOne,
                          a + p - 1 + (p - 1) * LDA, &intOne);
                    SWAP_(&k, a + k - 1, &LDA, a + p - 1, &LDA);
                    SWAP_(&kk, w + k - 1, &LDW, w + p - 1, &LDW);
                }
                if (kp != kk) {
                    a[kp - 1 + (k - 1) * LDA] = a[kk - 1 + (k - 1) * LDA];
                    intTmp = kp - k - 1;
                    COPY_(&intTmp, a + k + (kk - 1) * LDA, &intOne,
                          a + kp - 1 + k * LDA, &LDA);
                    intTmp = N - kp + 1;
                    COPY_(&intTmp, a + kp - 1 + (kk - 1) * LDA, &intOne,
                          a + kp - 1 + (kp - 1) * LDA, &intOne);
                    SWAP_(&kk, a + kk - 1, &LDA, a + kp - 1, &LDA);
                    SWAP_(&kk, w + kk - 1, &LDW, w + kp - 1, &LDW);
                }
                if (kStep == 1) {
                    intTmp = N - k + 1;
                    COPY_(&intTmp, w + k - 1 + (k - 1) * LDW, &intOne,
                          a + k - 1 + (k - 1) * LDA, &intOne);
                    if (k < N) {
                        if (ABS_(a[k - 1 + (k - 1) * LDA]) >= sfMin) {
                            r1 = T_ONE / a[k - 1 + (k - 1) * LDA];
                            intTmp = N - k;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&r1, a + k + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &r1, a + k + (k - 1) * LDA, &intOne);
#endif
                        } else if (a[k - 1 + (k - 1) * LDA] != ZERO) {
#pragma omp parallel for private(ii)
                            for (ii = k + 1; ii <= N; ii++) {
                                a[ii - 1 + (k - 1) * LDA] =
                                    a[ii - 1 + (k - 1) * LDA] /
                                    a[k - 1 + (k - 1) * LDA];
                            }
                        }
                    }
                } else {
                    if (k < N - 1) {
                        d21 = w[k + (k - 1) * LDW];
                        d11 = w[k + k * LDW] / d21;
                        d22 = w[k - 1 + (k - 1) * LDW] / d21;
                        t = T_ONE / (d11 * d22 - T_ONE);
#pragma omp parallel for private(j)
                        for (j = k + 2; j <= N; j++) {
                            a[j - 1 + (k - 1) * LDA] =
                                t * ((d11 * w[j - 1 + (k - 1) * LDW] -
                                      w[j - 1 + k * LDW]) /
                                     d21);
                            a[j - 1 + k * LDA] =
                                t * ((d22 * w[j - 1 + k * LDW] -
                                      w[j - 1 + (k - 1) * LDW]) /
                                     d21);
                        }
                    }
                    a[k - 1 + (k - 1) * LDA] = w[k - 1 + (k - 1) * LDW];
                    a[k + (k - 1) * LDA] = w[k + (k - 1) * LDW];
                    a[k + k * LDA] = w[k + k * LDW];
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
#pragma omp parallel for private(j, jj, jb, intTmp, Tmp2)
        for (j = k; j <= N; j += NB) {
            jb = MIN(NB, N - j + 1);
            for (jj = j; jj <= j + jb - 1; jj++) {
                intTmp = j + jb - jj;
                Tmp2 = k - 1;
                GEMV_("N", &intTmp, &Tmp2, &NEG_CONE, a + jj - 1, &LDA,
                      w + jj - 1, &LDW, &CONE, a + jj - 1 + (jj - 1) * LDA,
                      &intOne);
            }
            if (j + jb <= N) {
                intTmp = N - j - jb + 1;
                Tmp2 = k - 1;
                GEMM_("N", "T", &intTmp, &jb, &Tmp2, &NEG_CONE, a + j + jb - 1,
                      &LDA, w + j - 1, &LDW, &CONE,
                      a + j + jb - 1 + (j - 1) * LDA, &LDA);
            }
        }
        j = k - 1;
        while (j >= 1) {
            kStep = 1;
            jp1 = 1;
            jj = j;
            jp2 = ipiv[j - 1];
            if (jp2 < 0) {
                jp2 = -jp2;
                j--;
                jp1 = -ipiv[j - 1];
                kStep = 2;
            }
            j--;
            if (jp2 != jj && j >= 1) {
                SWAP_(&j, a + jp2 - 1, &LDA, a + jj - 1, &LDA);
            }
            jj = j + 1;
            if (jp1 != jj && kStep == 2) {
                SWAP_(&j, a + jp1 - 1, &LDA, a + jj - 1, &LDA);
            }
        }
        *kb = k - 1;
    }
    return;
}
