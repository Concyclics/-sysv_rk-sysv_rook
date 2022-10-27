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

    /*
     *
     *  -- LAPACK computational routine --
     *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
     *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG
     Ltd..--
     *
     *     .. Scalar Arguments ..
           CHARACTER          UPLO
           INTEGER            INFO, LDA, N
     *     ..
     *     .. Array Arguments ..
           INTEGER            IPIV( * )
           REAL               A( LDA, * ), E( * )
     *     ..
     *
     *  =====================================================================
     */

    const int LDA = *lda;
    const int N = *n;

    int upper, done;
    int i, iMax, j, jMax, iTemp, k, kk, kp, kStep, p, ii;
    dataAccu absAkk, Alpha, colMax, rowMax, sTemp, sfMin;
    dataType d11, d12, d21, d22, t, wk, wkm1, wkp1, z;

    Alpha = (ONE + sqrt(SEVTEN)) / EIGHT;

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
        xerbla_("SYTF2_ROOK", info, 7);
        return;
    }

    /*
     *    Initialize Alpha for use in choosing pivot block size.
     *
     *    Alpha = ( one+sqrt( sevten ) ) / eight
     *
     *    Compute machine safe minimum
     */

    sfMin = lamch_("S");
    if (upper) {
        /*
         *  Factorize A as U*D*U**T using the upper triangle of A
         *
         *  K is the main loop index, decreasing from N to 1 in steps of
         *  1 or 2
         */

        k = N;

        while (k > 0) {
            kStep = 1;
            p = k;

            /*
             *    Determine rows and columns to be interchanged and whether
             *    a 1-by-1 or 2-by-2 pivot block will be used
             */

            absAkk = ABS_(A[k - 1 + (k - 1) * (LDA)]);

            /*
             *    iMax is the row-index of the largest off-diagonal element in
             *    column K, and colMax is its absolute value.
             *    Determine both colMax and iMax.
             */

            if (k > 1) {
                intTmp = k - 1;
                iMax = I_AMAX(&intTmp, A + (k - 1) * (LDA), &intOne);
                colMax = ABS_(A[iMax - 1 + (k - 1) * (LDA)]);
            } else {
                colMax = ZERO;
            }

            if (MAX(absAkk, colMax) == ZERO) {
                //  Column K is zero or underflow: set INFO and continue
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
            } else {
                /*
                 *   Test for interchange
                 *
                 *   Equivalent to testing for (used to handle NaN and Inf)
                 *   absAkk.GE.Alpha*colMax
                 */

                if (absAkk >= Alpha * colMax) {
                    /*
                     *   no interchange,
                     *   use 1-by-1 pivot block
                     */
                    kp = k;
                } else {
                    done = false;
                    // Loop until pivot found
                    while (!done) {
                        /*
                         *   Begin pivot search loop body
                         *
                         *   jMax is the column-index of the largest
                         * off-diagonal element in row iMax, and rowMax is its
                         * absolute value. Determine both rowMax and jMax.
                         */

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

                        /*
                         *   Equivalent to testing for (used to handle NaN and
                         * Inf) ABS( A( iMax, iMax ) ).GE.Alpha*rowMax
                         */
                        if ((ABS_(A[iMax - 1 + (iMax - 1) * (LDA)]) >=
                             Alpha * rowMax)) {
                            /*
                             *    interchange rows and columns K and iMax,
                             *    use 1-by-1 pivot block
                             */
                            kp = iMax;
                            done = true;
                        }
                        /*
                         *    Equivalent to testing for rowMax .EQ. colMax,
                         *    used to handle NaN and Inf
                         */
                        else if (p == jMax || rowMax <= colMax) {
                            /*
                             *    interchange rows and columns K+1 and iMax,
                             *    use 2-by-2 pivot block
                             */
                            kp = iMax;
                            kStep = 2;
                            done = true;
                        } else {
                            // Pivot NOT found, set variables and repeat
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                        }
                        // End pivot search loop body
                    }
                }
                /*
                 *    Swap TWO rows and TWO columns
                 *    First swap
                 */

                if (kStep == 2 && p != k) {
                    /*
                     *    Interchange rows and column K and P in the leading
                     *    submatrix A(1:k,1:k) if we have a 2-by-2 pivot
                     */
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

                // Second swap
                kk = k - kStep + 1;
                if (kp != kk) {
                    /*
                     *    Interchange rows and columns KK and KP in the leading
                     *    submatrix A(1:k,1:k)
                     */
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

                /*
                 *    Update the leading submatrix A(1:k,1:k)
                 */

                if (kStep == 1) {
                    /*
                     *   1-by-1 pivot block D(k): column k now holds
                     *
                     *   W(k) = U(k)*D(k)
                     *
                     *  where U(k) is the k-th column of U.
                     */
                    if (k > 1) {
                        /*
                         *    Perform a rank-1 update of A(1:k-1,1:k-1) and
                         * store U(k) in column k.
                         */
                        if (ABS_(A[k - 1 + (k - 1) * (LDA)]) >= sfMin) {
                            /*
                             *  Perform a rank-1 update of A(1:k-1,1:k-1) as
                             *  A := A - U(k)*D(k)*U(k)**T
                             *     = A - W(k)*1/D(k)*W(k)**T
                             */
                            d11 = T_ONE / A[k - 1 + (k - 1) * (LDA)];
                            d11 = -d11;
                            intTmp = k - 1;
                            SYR_("U", &intTmp, &d11, A + (k - 1) * (LDA),
                                 &intOne, A, &LDA);
                            d11 = -d11;
                            /*
                             *  Store U(k) in column k
                             */
                            intTmp = k - 1;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + (k - 1) * LDA, &intOne);
#endif
                        } else {
                            /*
                             *  store L(k) in column k
                             */
                            d11 = A[k - 1 + (k - 1) * (LDA)];
                            #pragma omp parallel for private(ii)
                            for (ii = 1; ii < k; ii++) {
                                A[ii - 1 + (k - 1) * (LDA)] /= d11;
                            }
                            /*
                             *  Perform a rank-1 update of A(k+1:n,k+1:n) as
                             *  A := A - U(k)*D(k)*U(k)**T
                             *     = A - W(k)*(1/D(k))*W(k)**T
                             *     = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T
                             */
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
                            // Store U((k) and U((k-1)) in columns k and k-1 for
                            // row j
                            A[j - 1 + (k - 1) * (LDA)] = wk / d12;
                            A[j - 1 + (k - 2) * (LDA)] = wkm1 / d12;
                        }
                    }
                }
                // End column K is nonsingular
            }
            // Store details of the interchanges in array IPIV
            if (kStep == 1) {
                ipiv[k - 1] = kp;
            } else {
                ipiv[k - 1] = -p;
                ipiv[k - 2] = -kp;
            }
            // Descrease K and return to the start of the main loop
            k = k - kStep;
        }
    } else {
        /*
         *  Factorize A as L*D*L**T using the lower triangle of A
         *
         *  Initialize the unused last entry of the subdiagonal array E.
         */

        k = 1;

        while (k <= N) {
            kStep = 1;
            p = k;

            /*
             *  Determine rows and columns to be interchanged and whether
             *  a 1-by-1 or 2-by-2 pivot block will be used
             */

            absAkk = ABS_(A[k - 1 + (k - 1) * (LDA)]);

            /*
             *  iMax is the row-index of the largest off-diagonal element in
             *  column K, and colMax is its absolute value.
             *  Determine both colMax and iMax.
             */

            if (k < N) {
                intTmp = N - k;
                iMax = k + I_AMAX(&intTmp, A + k + (k - 1) * (LDA), &intOne);
                colMax = ABS_(A[iMax - 1 + (k - 1) * (LDA)]);
            } else {
                colMax = ZERO;
            }

            if (MAX(absAkk, colMax) == ZERO) {
                /*
                 *  Column K is zero or underflow: set INFO and continue
                 */
                if (*info == 0) {
                    *info = k;
                }
                kp = k;
            } else {
                /*
                 *  Test for interchange
                 *
                 *  Equivalent to testing for (used to handle NaN and Inf)
                 *  absAkk.GE.Alpha*colMax
                 */

                if (absAkk >= Alpha * colMax) {
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
                         *  jMax is the column-index of the largest off-diagonal
                         *  element in row iMax, and rowMax is its absolute
                         * value. Determine both rowMax and jMax.
                         */

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
                            // Pivot NOT found, set variables and repeat
                            p = iMax;
                            colMax = rowMax;
                            iMax = jMax;
                        }
                    }
                }
                /*
                 *  Swap TWO rows and TWO columns
                 *
                 *  First swap
                 */
                if (kStep == 2 && p != k) {
                    /*
                     *  Interchange rows and column K and P in the trailing
                     *  submatrix A(k:n,k:n) if we have a 2-by-2 pivot
                     */
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

                // Second swap
                kk = k + kStep - 1;
                if (kp != kk) {
                    /*
                     *  Interchange rows and columns KK and KP in the trailing
                     *  submatrix A(k:n,k:n)
                     */
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

                // update the trailing submatrix

                if (kStep == 1) {
                    if (k < N) {
                        // Perform a rank-1 update of A(k+1:n,k+1:n) and store
                        // L(k) in column k
                        if (ABS_(A[k - 1 + (k - 1) * (LDA)]) >= sfMin) {
                            d11 = T_ONE / A[k - 1 + (k - 1) * (LDA)];
                            d11 = -d11;
                            intTmp = N - k;
                            SYR_("L", &intTmp, &d11, A + k + (k - 1) * (LDA),
                                 &intOne, A + k + (k) * (LDA), &LDA);
                            d11 = -d11;
                            // Store L(k) in column k
                            intTmp = N - k;
#if defined(COMPLEX) || defined(COMPLEX16)
                            SCAL_(&intTmp, (void*)&d11, A + k + (k - 1) * LDA,
                                  &intOne);
#else
                            SCAL_(&intTmp, &d11, A + k + (k - 1) * LDA,
                                  &intOne);
#endif
                        } else {
                            // Store L(k) in column k
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
                            // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row
                            // J

                            wk = t * (d11 * A[j - 1 + (k - 1) * (LDA)] -
                                      A[j - 1 + (k) * (LDA)]);
                            wkp1 = t * (d22 * A[j - 1 + (k) * (LDA)] -
                                        A[j - 1 + (k - 1) * (LDA)]);

                            // Perform a rank-2 update of A(k+2:n,k+2:n)
                            #pragma omp parallel for private(i)
                            for (i = j; i <= N; i++) {
                                A[i - 1 + (j - 1) * (LDA)] -=
                                    wk * A[i - 1 + (k - 1) * (LDA)] / d21 +
                                    wkp1 * A[i - 1 + (k) * (LDA)] / d21;
                            }
                            // Store L(k) and L(k+1) in cols k and k+1 for row J

                            A[j - 1 + (k - 1) * (LDA)] = wk / d21;
                            A[j - 1 + (k) * (LDA)] = wkp1 / d21;
                        }
                    }
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
            // Increase K and return to the start of the main loop
            k += kStep;
            // End main loop
        }
    }
    return;
}
