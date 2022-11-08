/*******************************************************************************
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The realization of LASYF_RK.
 * Author: linhouzhong
 * Create: 2022-06-26
 *******************************************************************************/

#include "LASYF_RK.h"

void LASYF_RK(const char* uplo,
              const int* n,
              const int* nb,
              int* kb,
              dataType* a,
              const int* lda,
              dataType* e,
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
    int k, kw, kStep, p, iMax, jMax, iTemp, kk, kp, kkw, jb;
    int j, ii, jj, flag;
    dataAccu sfMin, absAKK, colMax, rowMax, sTemp;
    dataType d11, d12, d21, d22, r1, t;
    *info = 0;
    const dataAccu ALPHA = (ONE + sqrt(SEVTEN)) / EIGHT;
    char sla[] = {'S'};
    // Compute machine safe minimum
    sfMin = lamch_(sla);
    if (KmlLsame(*uplo, 'U')) {
        e[0] = ZERO;
        for (k = N; (NB >= N || k > N - NB + 1) && k > 0; k -= kStep) {
            kw = NB + k - N;
            kStep = 1;
            p = k;
            // Copy column k of a to column kw of w and update it
            COPY_(&k, a + (k - 1) * LDA, &intOne, w + (kw - 1) * LDW, &intOne);
            if (k < N) {
                intTmp = N - k;
                GEMV_("N", &k, &intTmp, &NEG_CONE, a + k * LDA, &LDA,
                      w + k - 1 + kw * LDW, &LDW, &CONE, w + (kw - 1) * LDW,
                      &intOne);
            }
            /*
            Determine rows and columns to be interchanged and
            whether*a 1 - by - 1 or 2 - by - 2 pivot block will be used
            */
            absAKK = ABS_(w[k - 1 + (kw - 1) * LDW]);
            /*
            IMAX is the row-index of the largest off-diagonal element in
            column k, and COLMAX is its absolute value.
            Determine both COLMAX and IMAX.
            */
            if (k > 1) {
                intTmp = k - 1;
                iMax = I_AMAX(&intTmp, w + (kw - 1) * LDW, &intOne);
                colMax = ABS_(w[iMax - 1 + (kw - 1) * LDW]);
            } else {
                colMax = ZERO;
            }
            if (MAX(absAKK, colMax) == ZERO) {
                // Column k is zero or underflow: set info and continue
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                COPY_(&k, w + (kw - 1) * LDW, &intOne, a + (k - 1) * LDA,
                      &intOne);
                // Set e( k ) to zero
                if (k > 1) {
                    e[k - 1] = ZERO;
                }
            } else {
                // Test for interchange

                /*
                Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
                (used to handle NaN and Inf)
                */
                if (absAKK >= ALPHA * colMax) {
                    // no interchange, use 1-by-1 pivot block
                    kp = k;
                } else {
                    flag = 1;
                    // Loop until pivot found
                    while (flag) {
                        // Copy column IMAX to column KW-1 of w and update
                        // it

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

                        /*
                        JMAX is the column-index of the largest
                        off-diagonal element in row IMAX, and ROWMAX is
                        its absolute value. Determine both ROWMAX and
                        JMAX.
                        */
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
                        /*
                        Equivalent to testing for ABS( w( IMAX, KW-1 )
                        ).GE.ALPHA*ROWMAX (used to handle NaN and Inf)
                        */
                        if (ABS_(w[iMax - 1 + (kw - 2) * LDW]) >=
                            ALPHA * rowMax) {
                            // interchange rows and columns K and
                            // IMAX,use 1-by-1 pivot block
                            kp = iMax;
                            // copy column KW-1 of w to column KW of w
                            COPY_(&k, w + (kw - 2) * LDW, &intOne,
                                  w + (kw - 1) * LDW, &intOne);
                            flag = 0;
                        } else if ((p == jMax) ||
                                   (rowMax <=
                                    colMax))  // Equivalent to testing
                                              // for
                                              // ROWMAX.EQ.COLMAX,(used
                                              // to handle NaN and Inf)
                        {
                            // interchange rows and columns K-1 and
                            // IMAX,use 2-by-2 pivot block
                            kp = iMax;
                            kStep = 2;
                            flag = 0;
                        } else {
                            // Pivot not found: set params and repeat
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                            // Copy updated JMAXth (next IMAXth) column
                            // to Kth of w
                            COPY_(&k, w + (kw - 2) * LDW, &intOne,
                                  w + (kw - 1) * LDW, &intOne);
                        }
                    }
                }

                kk = k - kStep + 1;
                // KKW is the column of w which corresponds to column KK of
                // a
                kkw = NB + kk - N;
                if ((kStep == 2) && (p != k)) {
                    // Copy non-updated column K to column P
                    intTmp = k - p;
                    COPY_(&intTmp, a + p + (k - 1) * LDA, &intOne,
                          a + p - 1 + p * LDA, &LDA);

                    COPY_(&p, a + (k - 1) * LDA, &intOne, a + (p - 1) * LDA,
                          &intOne);

                    // Interchange rows K and P in last N - K + 1 columns of
                    // a and last N-K+2 columns of w
                    intTmp = N - k + 1;
                    SWAP_(&intTmp, a + k - 1 + (k - 1) * LDA, &LDA,
                          a + p - 1 + (k - 1) * LDA, &LDA);

                    intTmp = N - kk + 1;
                    SWAP_(&intTmp, w + k - 1 + (kkw - 1) * LDW, &LDW,
                          w + p - 1 + (kkw - 1) * LDW, &LDW);
                }

                // Updated column kp is already stored in column KKW of w
                if (kp != kk) {
                    // Copy non-updated column KK to column kp

                    a[kp - 1 + (k - 1) * LDA] = a[kk - 1 + (k - 1) * LDA];

                    intTmp = k - 1 - kp;
                    COPY_(&intTmp, a + kp + (kk - 1) * LDA, &intOne,
                          a + kp - 1 + kp * LDA, &LDA);

                    COPY_(&kp, a + (kk - 1) * LDA, &intOne, a + (kp - 1) * LDA,
                          &intOne);

                    // Interchange rows KK and kp in last N-KK+1 columns
                    // of a and w
                    intTmp = N - kk + 1;
                    SWAP_(&intTmp, a + kk - 1 + (kk - 1) * LDA, &LDA,
                          a + kp - 1 + (kk - 1) * LDA, &LDA);
                    intTmp = N - kk + 1;
                    SWAP_(&intTmp, w + kk - 1 + (kkw - 1) * LDW, &LDW,
                          w + kp - 1 + (kkw - 1) * LDW, &LDW);
                }

                if (kStep == 1) {
                    /*
                    1-by-1 pivot block D(k): column KW of w now holds
                    w(k) = U(k)*D(k),where U(k) is the k-th column of U
                    */
                    // Store U(k) in column k of a
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
                            for (ii = 1; ii < k; ii++) {
                                a[ii - 1 + (k - 1) * LDA] /=
                                    a[k - 1 + (k - 1) * LDA];
                            }
                        }
                        // Store the superdiagonal element of D in array E
                        e[k - 1] = ZERO;
                    }
                } else {
                    /*
                    2-by-2 pivot block D(k): columns KW and KW-1 of w now
                    hold ( w(k-1) w(k) ) = ( U(k-1) U(k) )*D(k) where U(k)
                    and U(k-1) are the k-th and (k-1)-th columnsof U
                    */
                    if (k > 2) {
                        // Store U(k) and U(k-1) in columns k and k-1 of a
                        d12 = w[k - 2 + (kw - 1) * LDW];
                        d11 = w[k - 1 + (kw - 1) * LDW] / d12;
                        d22 = w[k - 2 + (kw - 2) * LDW] / d12;
                        t = T_ONE / (d11 * d22 - T_ONE);
#pragma omp parallel for private(j)
                        for (j = 1; j < k - 1; j++) {
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
                    /*
                    Copy diagonal elements of D(K) to a,
                    copy superdiagonal element of D(K) to E(K) and
                    ZERO out superdiagonal entry of a
                    */
                    a[k - 2 + (k - 2) * LDA] = w[k - 2 + (kw - 2) * LDW];
                    a[k - 2 + (k - 1) * LDA] = ZERO;
                    a[k - 1 + (k - 1) * LDA] = w[k - 1 + (kw - 1) * LDW];
                    e[k - 1] = w[k - 2 + (kw - 1) * LDW];
                    e[k - 2] = ZERO;
                }
            }
            // Store details of the interchanges in IPIV
            if (kStep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k - 2] = -kp;
            }
        }
        kw = NB + k - N;
/*
Update the upper triangle of A11 (= a(1:k,1:k)) as
A11 := A11 - U12*D*U12**T = A11 - U12*w**T
computing blocks of NB columns at a time
*/
#pragma omp parallel for private(j, jb, jj, intTmp, Tmp2)
        for (j = ((k - 1) / NB) * NB + 1; j >= 1; j -= NB) {
            jb = MIN(NB, k - j + 1);
            // Update the upper triangle of the diagonal block

            for (jj = j; jj < j + jb; jj++) {
                intTmp = jj - j + 1;
                Tmp2 = N - k;
                GEMV_("N", &intTmp, &Tmp2, &NEG_CONE, a + j - 1 + k * LDA, &LDA,
                      w + jj - 1 + kw * LDW, &LDW, &CONE,
                      a + j - 1 + (jj - 1) * LDA, &intOne);
            }
            // Update the rectangular superdiagonal block
            if (j >= 2) {
                intTmp = j - 1;
                Tmp2 = N - k;
                GEMM_("N", "T", &intTmp, &jb, &Tmp2, &NEG_CONE, a + k * LDA,
                      &LDA, w + j - 1 + kw * LDW, &LDW, &CONE,
                      a + (j - 1) * LDA, &LDA);
            }
        }
        // Set KB to the number of columns factorized
        *kb = N - k;
    } else {
        /*
        Factorize the leading columns of a using the lower triangle of a
        and working forwards, and compute the matrix w = L21*D for use in
        updating A22
        */

        // Initialize the unused last entry of the subdiagonal array E.
        e[N - 1] = ZERO;
        for (k = 1; (k < NB || NB >= N) && k <= N; k += kStep) {
            kStep = 1;
            p = k;
            // Copy column K of a to column K of w and update it
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
            /*
            Determine rows and columns to be interchanged and whether
            a 1-by-1 or 2-by-2 pivot block will be used
            */
            absAKK = ABS_(w[k - 1 + (k - 1) * LDW]);

            /*
            IMAX is the row-index of the largest off-diagonal element in
            column K, and COLMAX is its absolute value.
            */

            // Determine both COLMAX and IMAX.
            if (k < N) {
                intTmp = N - k;
                iMax = k + I_AMAX(&intTmp, w + k + (k - 1) * LDW, &intOne);
                colMax = ABS_(w[iMax - 1 + (k - 1) * LDW]);
            } else {
                colMax = ZERO;
            }
            if (MAX(absAKK, colMax) == ZERO) {
                // Column K is zero or underflow: set INFO and continue
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                intTmp = N - k + 1;
                COPY_(&intTmp, w + k - 1 + (k - 1) * LDW, &intOne,
                      a + k - 1 + (k - 1) * LDA, &intOne);
                // Set E( K ) to zero
                if (k < N) {
                    e[k - 1] = ZERO;
                }
            } else {
                /*
                Test for interchange
                Equivalent to testing for ABSAKK.GE.ALPHA*COLMAX
                (used to handle NaN and Inf)
                */
                if (absAKK >= ALPHA * colMax) {
                    // no interchange, use 1-by-1 pivot block
                    kp = k;
                } else {
                    flag = 1;
                    // Loop until pivot found
                    while (flag) {
                        // Copy column IMAX to column K+1 of w and update it

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
                        /*
                        JMAX is the column-index of the largest off-diagonal
                        element in row IMAX, and ROWMAX is its absolute
                        value.
                        */

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
                        /*
                        Equivalent to testing for ABS( w( IMAX, K+1 )
                        ).GE.ALPHA*ROWMAX (used to handle NaN and Inf)
                        */
                        if (ABS_(w[iMax - 1 + k * LDW]) >= ALPHA * rowMax) {
                            // interchange rows and columns K and IMAX,use
                            // 1-by-1 pivot block
                            kp = iMax;
                            // copy column K+1 of w to column K of w
                            intTmp = N - k + 1;
                            COPY_(&intTmp, w + k - 1 + k * LDW, &intOne,
                                  w + k - 1 + (k - 1) * LDW, &intOne);
                            flag = 0;
                        }
                        /*
                        Equivalent to testing for ROWMAX.EQ.COLMAX,
                        (used to handle NaN and Inf)
                        */
                        else if ((p == jMax) || (rowMax <= colMax)) {
                            // interchange rows and columns K+1 and IMAX,use
                            // 2-by-2 pivot block
                            kp = iMax;
                            kStep = 2;
                            flag = 0;
                        } else {
                            // Pivot not found: set params and repeat
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                            // Copy updated JMAXth (next IMAXth) column to
                            // Kth of w
                            intTmp = N - k + 1;
                            COPY_(&intTmp, w + k - 1 + k * LDW, &intOne,
                                  w + k - 1 + (k - 1) * LDW, &intOne);
                        }
                    }
                }
                kk = k + kStep - 1;
                if ((kStep == 2) && (p != k)) {
                    // Copy non-updated column K to column P
                    intTmp = p - k;
                    COPY_(&intTmp, a + k - 1 + (k - 1) * LDA, &intOne,
                          a + p - 1 + (k - 1) * LDA, &LDA);
                    intTmp = N - p + 1;
                    COPY_(&intTmp, a + p - 1 + (k - 1) * LDA, &intOne,
                          a + p - 1 + (p - 1) * LDA, &intOne);
                    // Interchange rows K and P in first K columns of a and
                    // first K + 1 columns of w
                    SWAP_(&k, a + k - 1, &LDA, a + p - 1, &LDA);
                    SWAP_(&kk, w + k - 1, &LDW, w + p - 1, &LDW);
                }
                // Updated column kp is already stored in column KK of w
                if (kp != kk) {
                    // Copy non-updated column KK to column KP
                    a[kp - 1 + (k - 1) * LDA] = a[kk - 1 + (k - 1) * LDA];
                    intTmp = kp - k - 1;
                    COPY_(&intTmp, a + k + (kk - 1) * LDA, &intOne,
                          a + kp - 1 + k * LDA, &LDA);
                    intTmp = N - kp + 1;
                    COPY_(&intTmp, a + kp - 1 + (kk - 1) * LDA, &intOne,
                          a + kp - 1 + (kp - 1) * LDA, &intOne);
                    // Interchange rows KK and KP in first KK columns of A
                    // and W
                    SWAP_(&kk, a + kk - 1, &LDA, a + kp - 1, &LDA);
                    SWAP_(&kk, w + kk - 1, &LDW, w + kp - 1, &LDW);
                }
                if (kStep == 1) {
                    /*
                    1-by-1 pivot block D(k): column k of w now holds
                    w(k) = L(k)*D(k)
                    where L(k) is the k-th column of L
                    Store L(k) in column k of a
                    */
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
#pragma omp parallel for schedule(static) private(ii)
                            for (ii = k; ii < N; ii++) {
                                a[ii + (k - 1) * LDA] /=
                                    a[k - 1 + (k - 1) * LDA];
                            }
                        }
                        e[k - 1] = ZERO;
                    }
                } else {
                    /*
                    2-by-2 pivot block D(k): columns k and k+1 of w now hold
                    ( w(k) w(k+1) ) = ( L(k) L(k+1) )*D(k)
                    where L(k) and L(k+1) are the k-th and (k+1)-th columns
                    of L
                    */
                    if (k < N - 1) {
                        // Store L(k) and L(k+1) in columns k and k+1 of a
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
                    /*
                    Copy diagonal elements of D(K) to a,
                    copy subdiagonal element of D(K) to E(K) and
                    ZERO out subdiagonal entry of a
                    */
                    a[k - 1 + (k - 1) * LDA] = w[k - 1 + (k - 1) * LDW];
                    a[k + (k - 1) * LDA] = ZERO;
                    a[k + k * LDA] = w[k + k * LDW];
                    e[k - 1] = w[k + (k - 1) * LDW];
                    e[k] = ZERO;
                }
                // End column K is nonsingular
            }
            // Store details of the interchanges in IPIV
            if (kStep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k] = -kp;
            }
        }
/*
Update the lower triangle of A22 (= a(k:N,k:N)) as
A22 := A22 - L21*D*L21**T = A22 - L21*w**T
computing blocks of NB columns at a time
*/
#pragma omp parallel for private(j, jj, jb, intTmp, Tmp2)
        for (j = k; j <= N; j += NB) {
            jb = MIN(NB, N - j + 1);
            // Update the lower triangle of the diagonal block
            for (jj = j; jj <= j + jb - 1; jj++) {
                intTmp = j + jb - jj;
                Tmp2 = k - 1;
                GEMV_("N", &intTmp, &Tmp2, &NEG_CONE, a + jj - 1, &LDA,
                      w + jj - 1, &LDW, &CONE, a + jj - 1 + (jj - 1) * LDA,
                      &intOne);
            }
            // Update the rectangular subdiagonal block
            if (j + jb <= N) {
                intTmp = N - j - jb + 1;
                Tmp2 = k - 1;
                GEMM_("N", "T", &intTmp, &jb, &Tmp2, &NEG_CONE, a + j + jb - 1,
                      &LDA, w + j - 1, &LDW, &CONE,
                      a + j + jb - 1 + (j - 1) * LDA, &LDA);
            }
        }

        // Set KB to the number of columns factorized
        *kb = k - 1;
    }
    return;
}
