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

            /*
             *    Determine rows and columns to be interchanged and whether
             *    a 1-by-1 or 2-by-2 pivot block will be used
             */

            absakk = ABS_(A[k - 1 + (k - 1) * LDA]);

            /*
             *    IMAX is the row-index of the largest off-diagonal element in
             *    column K, and COLMAX is its absolute value.
             *    Determine both COLMAX and IMAX.
             */

            if (k > 1) {
                intTmp = k - 1;
                imax = I_AMAX(&intTmp, A + (k - 1) * LDA, &intOne);
                colmax = ABS_(A[imax - 1 + (k - 1) * LDA]);
            } else {
                colmax = ZERO;
            }

            if (MAX(absakk, colmax) == ZERO) {
                // printf("absakk = %f, colmax = %f MAX(absakk, colmax) = %f
                // MAX(absakk, colmax) == ZERO = %d\n", absakk, colmax,
                // MAX(absakk, colmax), (MAX(absakk, colmax) == ZERO));
                //  Column K is zero or underflow: set INFO and continue
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                if (k > 1) {
                    // Set E( K ) to zero
                    E[k - 1] = ZERO;
                }
            } else {
                /*
                 *   Test for interchange
                 *
                 *   Equivalent to testing for (used to handle NaN and Inf)
                 *   ABSAKK.GE.ALPHA*COLMAX
                 */

                if (absakk >= ALPHA * colmax) {
                    kp = k;
                } else {
                    done = false;
                    // Loop until pivot found

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

                        /*
                         *   Equivalent to testing for (used to handle
                         * NaN and Inf) ABS( A( IMAX, IMAX )
                         * ).GE.ALPHA*ROWMAX
                         */
                        if ((ABS_(A[imax - 1 + (imax - 1) * LDA]) >=
                             ALPHA * rowmax)) {
                            /*
                             *    interchange rows and columns K and
                             * IMAX, use 1-by-1 pivot block
                             */
                            kp = imax;
                            done = true;
                        }
                        /*
                         *    Equivalent to testing for ROWMAX .EQ.
                         * COLMAX, used to handle NaN and Inf
                         */
                        else if (p == jmax || rowmax <= colmax) {
                            /*
                             *    interchange rows and columns K+1 and
                             * IMAX, use 2-by-2 pivot block
                             */
                            kp = imax;
                            done = true;
                            kstep = 2;
                        } else {
                            // Pivot NOT found, set variables and repeat
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                        }
                        // End pivot search loop body
                    }
                }
                /*
                 *    Swap TWO rows and TWO columns
                 *    First swap
                 */
                if (kstep == 2 && p != k) {
                    /*
                     *    Interchange rows and column K and P in the leading
                     *    submatrix A(1:k,1:k) if we have a 2-by-2 pivot
                     */
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
                    /*
                     *    Convert upper triangle of A into U form by applying
                     *    the interchanges in columns k+1:N.
                     */
                    if (k < N) {
                        intTmp = N - k;
                        SWAP_(&intTmp, A + k - 1 + k * LDA, &LDA,
                              A + p - 1 + k * LDA, &LDA);
                    }
                }

                // Second swap
                kk = k - kstep + 1;
                if (kp != kk) {
                    // printf("kp = %d, kk = %d\n", kp, kk);
                    /*
                     *    Interchange rows and columns KK and KP in the leading
                     *    submatrix A(1:k,1:k)
                     */
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
                    /*
                     *    Convert upper triangle of A into U form by applying
                     *    the interchanges in columns k+1:N.
                     */
                    if (k < N) {
                        intTmp = N - k;
                        SWAP_(&intTmp, A + kk - 1 + k * LDA, &LDA,
                              A + kp - 1 + k * LDA, &LDA);
                    }
                }
                if (kstep == 1) {
                    // printf("kstep = %d, p = %d, k = %d\n", kstep, p, k);
                    /*
                     *   1-by-1 pivot block Dk: column k now holds
                     *
                     *   Wk = Uk*Dk
                     *
                     *  where Uk is the k-th column of U.
                     */
                    if (k > 1) {
                        /*
                         *    Perform a rank-1 update of A(1:k-1,1:k-1) and
                         * store Uk in column k.
                         */
                        if (ABS_(A[k - 1 + (k - 1) * LDA]) >= sfmin) {
                            /*
                             *  Perform a rank-1 update of A(1:k-1,1:k-1) as
                             *  A := A - Uk*Dk*Uk**T
                             *     = A - Wk*1/Dk*Wk**T
                             */
                            d11 = T_ONE / A[k - 1 + (k - 1) * LDA];
                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * LDA, &intOne,
                                 A, &LDA);
                            d11 = -d11;
                            // Store Uk in column k
                            intTmp = k - 1;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + (k - 1) * LDA, &intOne);
#endif
                        } else {
                            // store Lk in column k
                            d11 = A[k - 1 + (k - 1) * LDA];

#pragma omp parallel for private(ii)
                            for (ii = 1; ii < k; ii++) {
                                A[ii - 1 + (k - 1) * LDA] /= d11;
                            }
                            /*
                             *  Perform a rank-1 update of A(k+1:n,k+1:n) as
                             *  A := A - Uk*Dk*Uk**T
                             *     = A - Wk*(1/Dk)*Wk**T
                             *     = A - (Wk/Dk)*(Dk)*(Wk/D(K))**T
                             */
                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * LDA, &intOne,
                                 A, &LDA);
                            d11 = -d11;
                        }
                        // Store the superdiagonal element of D in array E
                        E[k - 1] = ZERO;
                    }
                } else {
                    /*
                     *  2-by-2 pivot block Dk: columns k and k-1 now hold
                     *
                     *  ( W(k-1) Wk ) = ( U(k-1) Uk )*Dk
                     *
                     *  where Uk and U(k-1) are the k-th and (k-1)-th columns
                     *  of U
                     *
                     *  Perform a rank-2 update of A(1:k-2,1:k-2) as
                     *
                     *  A := A - ( U(k-1) Uk )*Dk*( U(k-1) Uk )**T
                     *     = A - ( ( A(k-1)Ak )*inv(Dk) ) * ( A(k-1)Ak
                     * )**T
                     *
                     *  and store Lk and L(k+1) in columns k and k+1
                     */
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
                            // Store U(k and U((k-1)) in columns k and k-1 for
                            // row j
                            A[j - 1 + (k - 1) * LDA] = wk / d12;
                            A[j - 1 + (k - 2) * LDA] = wkm1 / d12;
                        }
                    }
                    /*
                     *  Copy superdiagonal elements of D(K) to E(K) and
                     *  ZERO out superdiagonal entry of A
                     */
                    E[k - 1] = A[k - 2 + (k - 1) * LDA];
                    E[k - 2] = ZERO;
                    A[k - 2 + (k - 1) * LDA] = ZERO;
                }
                // End column K is nonsingular
            }
            // Store details of the interchanges in array IPIV
            if (kstep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k - 2] = -kp;
            }
            // Descrease K and return to the start of the main loop
            k = k - kstep;
        }
    } else {
        /*
         *  Factorize A as L*D*L**T using the lower triangle of A
         *
         *  Initialize the unused last entry of the subdiagonal array E.
         */

        E[N - 1] = ZERO;

        // K is the main loop index, increasing from 1 to N in steps of 1 or 2

        k = 1;

        while (k <= N) {
            kstep = 1;
            p = k;

            /*
             *  Determine rows and columns to be interchanged and whether
             *  a 1-by-1 or 2-by-2 pivot block will be used
             */

            absakk = ABS_(A[k - 1 + (k - 1) * LDA]);

            /*
             *  IMAX is the row-index of the largest off-diagonal element in
             *  column K, and COLMAX is its absolute value.
             *  Determine both COLMAX and IMAX.
             */

            if (k < N) {
                intTmp = N - k;
                imax = k + I_AMAX(&intTmp, A + k + (k - 1) * LDA, &intOne);
                colmax = ABS_(A[imax - 1 + (k - 1) * LDA]);
            } else {
                colmax = ZERO;
            }

            if (MAX(absakk, colmax) == ZERO) {
                /*
                 *  Column K is zero or underflow: set INFO and continue
                 */
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
                // Set E(K) to zero if K < N
                if (k < N) {
                    E[k - 1] = ZERO;
                }
            } else {
                /*
                 *  Test for interchange
                 *
                 *  Equivalent to testing for (used to handle NaN and Inf)
                 *  ABSAKK.GE.ALPHA*COLMAX
                 */

                if (absakk >= ALPHA * colmax) {
                    /*
                     *  No interchange, use 1-by-1 pivot block
                     */
                    kp = k;
                } else {
                    done = false;
                    while (!done) {
                        /*
                         *  Begin pivot search loop body
                         *
                         *  JMAX is the column-index of the largest off-diagonal
                         *  element in row IMAX, and ROWMAX is its absolute
                         * value. Determine both ROWMAX and JMAX.
                         */

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
                        /*
                         *  Equivalent to testing for (used to handle NaN and
                         * Inf) ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX
                         */
                        if (ABS_(A[imax - 1 + (imax - 1) * LDA]) >=
                            ALPHA * rowmax) {
                            /*
                             *  interchange rows and columns K and IMAX,
                             *  use 1-by-1 pivot block
                             */
                            kp = imax;
                            done = true;
                        } else if (p == jmax || rowmax <= colmax) {
                            /*
                             *  Equivalent to testing for ROWMAX .EQ. COLMAX,
                             *  used to handle NaN and Inf
                             *
                             *  interchange rows and columns K+1 and IMAX,
                             *  use 2-by-2 pivot block
                             */
                            kp = imax;
                            kstep = 2;
                            done = true;
                        } else {
                            // Pivot NOT found, set variables and repeat
                            p = imax;
                            colmax = rowmax;
                            imax = jmax;
                        }
                    }
                }
                /*
                 *  Swap TWO rows and TWO columns
                 *
                 *  First swap
                 */
                if (kstep == 2 && p != k) {
                    /*
                     *  Interchange rows and column K and P in the trailing
                     *  submatrix A(k:n,k:n) if we have a 2-by-2 pivot
                     */
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
                    /*
                     *  Convert lower triangle of A into L form by applying
                     *  the interchanges in columns 1:k-1.
                     */
                    if (k > 1) {
                        intTmp = k - 1;
                        SWAP_(&intTmp, A + k - 1, &LDA, A + p - 1, &LDA);
                    }
                }

                // Second swap
                kk = k + kstep - 1;
                if (kp != kk) {
                    /*
                     *  Interchange rows and columns KK and KP in the trailing
                     *  submatrix A(k:n,k:n)
                     */
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
                    /*
                     *  Convert lower triangle of A into L form by applying
                     *  the interchanges in columns 1:k-1.
                     */
                    if (k > 1) {
                        intTmp = k - 1;
                        SWAP_(&intTmp, A + kk - 1, &LDA, A + kp - 1, &LDA);
                    }
                }

                // update the trailing submatrix

                if (kstep == 1) {
                    /*
                     *  1-by-1 pivot block Dk: column k now holds
                     *
                     *  Wk = Lk*Dk
                     *
                     *  where Lk is the k-th column of L
                     */
                    if (k < N) {
                        // Perform a rank-1 update of A(k+1:n,k+1:n) and store
                        // Lk in column k
                        if (ABS_(A[k - 1 + (k - 1) * LDA]) >= sfmin) {
                            /*
                             *  Perform a rank-1 update of A(k+1:n,k+1:n) as
                             *  A := A - Lk*Dk*Lk**T
                             *     = A - Wk*(1/Dk)*Wk**T
                             */
                            d11 = T_ONE / A[k - 1 + (k - 1) * LDA];
                            intTmp = N - k;
                            d11 = -d11;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * LDA,
                                 &intOne, A + k + k * LDA, &LDA);
                            d11 = -d11;
                            // Store Lk in column k
                            intTmp = N - k;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + k + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + k + (k - 1) * LDA,
                                  &intOne);
#endif
                        } else {
                            // Store Lk in column k
                            d11 = A[k - 1 + (k - 1) * LDA];

#pragma omp parallel for private(ii)
                            for (ii = k + 1; ii <= N; ii++) {
                                A[ii - 1 + (k - 1) * LDA] /= d11;
                            }
                            /*
                             *  Perform a rank-1 update of A(k+1:n,k+1:n) as
                             *  A := A - Lk*Dk*Lk**T
                             *     = A - Wk*(1/Dk)*Wk**T
                             *     = A - (Wk/Dk)*(Dk)*(Wk/D(K))**T
                             */
                            d11 = -d11;
                            intTmp = N - k;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * LDA,
                                 &intOne, A + k + k * LDA, &LDA);
                            d11 = -d11;
                        }
                        // Store the subdiagonal element of D in array E
                        E[k - 1] = ZERO;
                    }
                } else {
                    /*
                     *  2-by-2 pivot block Dk: columns k and k+1 now hold
                     *
                     *  ( Wk W(k+1) ) = ( Lk L(k+1) )*Dk
                     *
                     *  where Lk and L(k+1) are the k-th and (k+1)-th columns
                     *  of L
                     *
                     *
                     *  Perform a rank-2 update of A(k+2:n,k+2:n) as
                     *
                     *  A := A - ( Lk L(k+1) ) * Dk * ( Lk L(k+1) )**T
                     *     = A - ( ( AkA(k+1) )*inv(Dk ) * ( AkA(k+1) )**T
                     *
                     *  and store Lk and L(k+1) in columns k and k+1
                     */
                    if (k < N - 1) {
                        d21 = A[k + (k - 1) * LDA];
                        d11 = A[k + k * LDA] / d21;
                        d22 = A[k - 1 + (k - 1) * LDA] / d21;
                        t = T_ONE / (d11 * d22 - T_ONE);

                        for (j = k + 2; j <= N; ++j) {
                            // Compute  D21 * ( WkW(k+1) ) * inv(Dk) for row
                            // J

                            wk = t * (d11 * A[j - 1 + (k - 1) * LDA] -
                                      A[j - 1 + k * LDA]);
                            wkp1 = t * (d22 * A[j - 1 + k * LDA] -
                                        A[j - 1 + (k - 1) * LDA]);

                            // Perform a rank-2 update of A(k+2:n,k+2:n)
#pragma omp parallel for private(i)
                            for (i = j; i <= N; ++i) {
                                A[i - 1 + (j - 1) * LDA] -=
                                    wk * A[i - 1 + (k - 1) * LDA] / d21 +
                                    wkp1 * A[i - 1 + k * LDA] / d21;
                            }
                            // Store Lk and L(k+1) in cols k and k+1 for row J

                            A[j - 1 + (k - 1) * LDA] = wk / d21;
                            A[j - 1 + k * LDA] = wkp1 / d21;
                        }
                    }
                    /*
                     *  Copy subdiagonal elements of D(K) to E(K) and
                     *  ZERO out subdiagonal entry of A
                     */
                    E[k - 1] = A[k + (k - 1) * LDA];
                    E[k] = ZERO;
                    A[k + (k - 1) * LDA] = ZERO;
                }
                // End column K is nonsingular
            }
            // Store details of the interchanges in IPIV
            if (kstep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k] = -kp;
            }
            // Increase K and return to the start of the main loop
            k += kstep;
            // End main loop
        }
    }
    return;
}
