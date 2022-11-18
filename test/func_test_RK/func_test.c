/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of _FuncTest
 * Author: linhouzhong
 * Create: 2022-06-26
 */

#include "func_test.h"

double calc_rcond(char uplo,
                  int n,
                  funcTdataType* A,
                  int lda,
                int* ipiv){
    funcTdataAccu anorm;
    funcTdataAccu rcond;
    int info;
    funcTdataAccu work[2 * n];
    funcTdataType W1[2 * n];
    int iwork[n];
    anorm = _lansy_("I", &uplo, &n, A, &lda, work);
    #if defined(SINGLE)
        ssycon_(&uplo, &n, A, &lda, ipiv, &anorm, &rcond, W1, iwork, &info);
    #elif defined(DOUBLE)
        dsycon_(&uplo, &n, A, &lda, ipiv, &anorm, &rcond, W1, iwork, &info);
    #elif defined(COMPLEX)
        csycon_(&uplo, &n, A, &lda, ipiv, &anorm, &rcond, W1, &info);
    #elif defined(COMPLEX16)
        zsycon_(&uplo, &n, A, &lda, ipiv, &anorm, &rcond, W1, &info);
    #endif
    return rcond;
}

double _TestWithSeed(int rand_seed,
                     int N,
                     int NRHS,
                     char uplo,
                     funcTdataType* A,
                     funcTdataType* A1,
                     funcTdataType* E,
                     funcTdataType* B1,
                     funcTdataType* B3,
                     funcTdataType* W,
                     funcTdataAccu* W1,
                     int* ipiv) {
    srand(rand_seed);

    int size = sizeof(funcTdataType);

    int info;
    funcTdataAccu RESID;
    int lwork = N * 64;

    // printf("  UPLO = %c\n", uplo);
    _CreatMatrix(uplo, N, NRHS, A, B1);
    // copy A to A1 and A2.
    for (int j = 0; j < N * N; ++j) {
        A1[j] = A[j];
    }
    // copy B to B1 and B2.
    for (int j = 0; j < N * NRHS; ++j) {
        B3[j] = B1[j];
    }

    // calculate result by lapack code and our code
    SYSV_RK(&uplo, &N, &NRHS, A1, &N, E, ipiv, B1, &N, W, &lwork,
            &info);  // lapack code
    _pot02_(&uplo, &N, &NRHS, A, &N, B1, &N, B3, &N, W1,
            &RESID);  // lapack code
    return RESID;
}

double _TestWithSeedGet04(int rand_seed,
                          int N,
                          int NRHS,
                          char uplo,
                          funcTdataType* A,
                          funcTdataType* A1,
                          funcTdataType* E,
                          funcTdataType* XACT,
                          funcTdataType* B,
                          funcTdataType* W,
                          funcTdataAccu* W1,
                          int* ipiv) {
    srand(rand_seed);

    int size = sizeof(funcTdataType);

    int info;
    funcTdataAccu RESID, RCOND;
    int lwork = N * 64;

    // printf("  UPLO = %c\n", uplo);
    _CreatMatrix(uplo, N, NRHS, A, XACT);

    funcTdataType* Af = (funcTdataType*)malloc(N * N * size);
    for (int j = 0; j < N * N; ++j) {
        Af[j] = A[j];
    }

    _sytrf_(&uplo, &N, Af, &N, ipiv, W, &N, &info);
    printf("info = %d\n", info);
    RCOND = calc_rcond(uplo, N, Af, N, ipiv);
    free(Af);

    if (uplo == 'L') {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < i; ++j) {
                A[i * N + j] = A[j * N + i];
            }
        }
    } else {
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                A[i * N + j] = A[j * N + i];
            }
        }
    }

    // B = A*X
    funcTdataType alpha = 1.0;
    funcTdataType beta = 0.0;

    // copy A to A1.
    for (int j = 0; j < N * N; ++j) {
        A1[j] = A[j];
    }

     _gemm_("N", "N", &N, &NRHS, &N, &alpha, A, &N, XACT, &N, &beta, B, &N);

    SYSV_RK(&uplo, &N, &NRHS, A1, &N, E, ipiv, B, &N, W, &lwork,
            &info);  // lapack code

    _get04_(&N, &NRHS, B, &N, XACT, &N, &RCOND, &RESID);  // lapack code

    return RESID;
}

int cover() {
    int N, NRHS, NB;
    char UPLO;
    N = 20;
    NRHS = 20;
    NB = 16;
    UPLO = 'U';
    funcTdataType *A, *A1, *E, *B1, *B3, *W;
    funcTdataAccu* W1;
    int* ipiv;
    char sla[] = {'S'};
    dataAccu sfMin = lamch_(sla);
    int size = sizeof(funcTdataType);
    if (((A = (funcTdataType*)malloc((long int)size * N * N)) == NULL) ||
        ((A1 = (funcTdataType*)malloc((long int)size * N * N)) == NULL) ||
        ((E = (funcTdataType*)malloc((long int)size * N)) == NULL) ||
        ((B1 = (funcTdataType*)malloc((long int)size * N * NRHS)) == NULL) ||
        ((B3 = (funcTdataType*)malloc((long int)size * N * NRHS)) == NULL) ||
        ((W = (funcTdataType*)malloc((long int)size * N * N)) == NULL) ||
        ((W1 = (funcTdataAccu*)malloc((long int)sizeof(funcTdataAccu) * N *
                                      N)) == NULL) ||
        ((ipiv = (int*)malloc((long int)sizeof(int) * N)) == NULL)) {
        printf("malloc error!\n");
        return -1;
    }

    int info;
    // A[all] = sfmin/2;
    UPLO = 'U';
    for (int i = 0; i < N * N; i++) {
        A[i] = sfMin / 2;
    }
    LASYF_RK(&UPLO, &N, &NB, &N, A, &N, E, ipiv, W, &N, &info);
    for (int i = 0; i < N * N; i++) {
        A[i] = sfMin / 2;
    }
    SYTF2_RK(&UPLO, &N, A, &N, E, ipiv, &info);
    UPLO = 'L';
    for (int i = 0; i < N * N; i++) {
        A[i] = sfMin / 2;
    }
    LASYF_RK(&UPLO, &N, &NB, &N, A, &N, E, ipiv, W, &N, &info);
    for (int i = 0; i < N * N; i++) {
        A[i] = sfMin / 2;
    }
    SYTF2_RK(&UPLO, &N, A, &N, E, ipiv, &info);
    // A[all]=0
    UPLO = 'U';
    for (int i = 0; i < N * N; i++) {
        A[i] = 0;
    }
    SYSV_RK(&UPLO, &N, &NRHS, A, &N, E, ipiv, B1, &N, W, &N, &info);
    UPLO = 'L';
    for (int i = 0; i < N * N; i++) {
        A[i] = 0;
    }
    SYSV_RK(&UPLO, &N, &NRHS, A, &N, E, ipiv, B1, &N, W, &N, &info);
    // N=NRHS=0
    for (int i = 0; i < N * N; i++) {
        A[i] = rand();
    }
    int lda = 1, ldb = 1, ldw = 1;
    int zero_0 = 0, neg_1 = -1;
    UPLO = 'U';
    SYSV_RK(&UPLO, &zero_0, &zero_0, A, &lda, E, ipiv, B1, &ldb, W, &ldw,
            &info);
    UPLO = 'L';
    SYSV_RK(&UPLO, &zero_0, &zero_0, A, &lda, E, ipiv, B1, &ldb, W, &ldw,
            &info);
    // lwork<nb
    N = 20;
    SYTRF_RK(&UPLO, &N, A, &N, E, ipiv, W, &N, &info);
    // lwork=-1
    SYSV_RK(&UPLO, &N, &NRHS, A, &lda, E, ipiv, B1, &ldb, W, &neg_1, &info);
    // SYTRF parameter error
    char non_uplo = 'a';
    SYTRF_RK(&non_uplo, &N, A, &N, E, ipiv, W, &N, &info);
    SYTRF_RK(&UPLO, &neg_1, A, &N, E, ipiv, W, &N, &info);
    SYTRF_RK(&UPLO, &N, A, &neg_1, E, ipiv, W, &N, &info);
    SYTRF_RK(&UPLO, &N, A, &N, E, ipiv, W, &zero_0, &info);
    UPLO = 'U';
    SYTRF_RK(&UPLO, &zero_0, A, &lda, E, ipiv, W, &ldw, &info);
    UPLO = 'L';
    SYTRF_RK(&UPLO, &zero_0, A, &lda, E, ipiv, W, &ldw, &info);
    // SYTF2 parameter error
    SYTF2_RK(&non_uplo, &N, A, &N, E, ipiv, &info);
    SYTF2_RK(&UPLO, &neg_1, A, &N, E, ipiv, &info);
    SYTF2_RK(&UPLO, &N, A, &neg_1, E, ipiv, &info);
    // SYTRS parameter error
    SYTRS_3(&non_uplo, &N, &NRHS, A, &N, E, ipiv, B1, &N, &info);
    SYTRS_3(&UPLO, &neg_1, &NRHS, A, &N, E, ipiv, B1, &N, &info);
    SYTRS_3(&UPLO, &N, &neg_1, A, &N, E, ipiv, B1, &N, &info);
    SYTRS_3(&UPLO, &N, &NRHS, A, &neg_1, E, ipiv, B1, &N, &info);
    SYTRS_3(&UPLO, &N, &NRHS, A, &N, E, ipiv, B1, &neg_1, &info);

    free(A);
    free(A1);
    free(E);
    free(B1);
    free(B3);
    free(W);
    free(W1);
    free(ipiv);
    return 0;
}

void _FuncTest(int times,
               char scale,
               funcTdataType* A,
               funcTdataType* A1,
               funcTdataType* E,
               funcTdataType* B1,
               funcTdataType* B2,
               funcTdataType* B3,
               funcTdataType* B4,
               funcTdataType* W,
               funcTdataAccu* W1,
               int* ipiv) {
    srand((unsigned)time(NULL));

    int size = sizeof(funcTdataType);

    for (int i = 0; i < times; ++i) {
        int N, NRHS;
        int info1, info2;
        funcTdataAccu RESID1, RESID2;
        char uplo = ((rand() % 2) == 1) ? 'L' : 'U';
        printf("Test %d\n", i + 1);
        if (scale == 's' || scale == 'S') {
            N = (rand() % 991) + 10;  //[10,1000]
        } else if (scale == 'm' || scale == 'M') {
            N = (rand() % 9001) + 1000;  //[1000,10000]
        } else if (scale == 'b' || scale == 'B') {
            N = (rand() % 30001) + 10000;  //[10000,40000]
        }
        int lwork = N * 64;
        NRHS = rand() % N + 1;  // NRHS=[1,N]

        printf("  UPLO = %c\n", uplo);
        _CreatMatrix(uplo, N, NRHS, A, B1);
        // copy A to A1 and A2.
        for (int j = 0; j < N * N; ++j) {
            A1[j] = A[j];
        }
        // copy B to B1 and B2.
        for (int j = 0; j < N * NRHS; ++j) {
            funcTdataType temp = B1[j];
            B2[j] = temp;
            B3[j] = temp;
            B4[j] = temp;
        }

        // calculate result by lapack code and our code
        _sysv_rk_(&uplo, &N, &NRHS, A1, &N, E, ipiv, B1, &N, W, &lwork,
                  &info1);  // lapack code
        _pot02_(&uplo, &N, &NRHS, A, &N, B1, &N, B3, &N, W1,
                &RESID1);  // lapack code
        // copy A to A1 and A2.
        for (int j = 0; j < N * N; ++j) {
            A1[j] = A[j];
        }
        SYSV_RK(&uplo, &N, &NRHS, A1, &N, E, ipiv, B2, &N, W, &lwork,
                &info2);  // our code
        _pot02_(&uplo, &N, &NRHS, A, &N, B2, &N, B4, &N, W1,
                &RESID2);  // our code

        printf("  resid1 = %f,resid2 = %f\n", RESID1, RESID2);

        if (RESID1 >= RESID2) {
            printf("  N = %d, NRHS = %d,  succeed!\n\n", N, NRHS);
        } else {
            printf("  N = %d, NRHS = %d,  NO succeed...\n\n", N, NRHS);
        }
    }
}

void nullptr_test(int N, int M) {
    srand((unsigned)time(NULL));
    char uplo = ((rand() % 2) == 1) ? 'L' : 'U';
    funcTdataType *A, *B, *E, *W;
    int *ipiv, info;
    struct timespec start, end;
    double difTime, result;
    float opTRF, opTRS;
    int lin0 = 0;
    if (((A = (funcTdataType*)malloc(sizeof(funcTdataType) * N * N)) == NULL) ||
        ((E = (funcTdataType*)malloc(sizeof(funcTdataType) * N)) == NULL) ||
        ((B = (funcTdataType*)malloc(sizeof(funcTdataType) * N * M)) == NULL) ||
        ((W = (funcTdataType*)malloc(sizeof(funcTdataType) * N * N)) == NULL) ||
        ((ipiv = (int*)malloc(sizeof(int) * N)) == NULL)) {
        printf("malloc error!\n");
        return;
    }
    _CreatMatrix(uplo, N, M, A, B);
    // SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &N, W, &N, &info);	// our
    //  code
    int NRHS = N;
    // test nullptr
    printf("Test nullptr\n");
    printf("test nullptr for param 1, uplo\n");
    SYSV_RK(NULL, &N, &NRHS, A, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 2, N\n");
    SYSV_RK(&uplo, NULL, &NRHS, A, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 3, NRHS\n");
    SYSV_RK(&uplo, &N, NULL, A, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 4, A\n");
    SYSV_RK(&uplo, &N, &NRHS, NULL, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 5, lda\n");
    SYSV_RK(&uplo, &N, &NRHS, A, NULL, E, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 6, E\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, NULL, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 7, ipiv\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, NULL, B, &N, W, &N, &info);
    printf("test nullptr for param 8, B\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, NULL, &N, W, &N, &info);
    printf("test nullptr for param 9, ldb\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, NULL, W, &N, &info);
    printf("test nullptr for param 10, W\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &N, NULL, &N, &info);
    printf("test nullptr for param 11, ldwork\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &N, W, NULL, &info);
    printf("test nullptr for param 12, info\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &N, W, &N, NULL);

    free(A);
    free(B);
    free(E);
    free(W);
    free(ipiv);

    printf("test nullptr end!\n");
}

void exception_test(int N, int M) {
    srand((unsigned)time(NULL));
    char uplo = ((rand() % 2) == 1) ? 'L' : 'U';
    funcTdataType *A, *B, *E, *W;
    int *ipiv, info;
    if (((A = (funcTdataType*)malloc(sizeof(funcTdataType) * N * N)) == NULL) ||
        ((E = (funcTdataType*)malloc(sizeof(funcTdataType) * N)) == NULL) ||
        ((B = (funcTdataType*)malloc(sizeof(funcTdataType) * N * M)) == NULL) ||
        ((W = (funcTdataType*)malloc(sizeof(funcTdataType) * N * N)) == NULL) ||
        ((ipiv = (int*)malloc(sizeof(int) * N)) == NULL)) {
        printf("malloc error!\n");
        return;
    }
    _CreatMatrix(uplo, N, M, A, B);
    // SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &N, W, &N, &info);	// our
    //  code
    int NRHS = N;
    int neg = -2;
    char wronguplo = 'a';
    printf("test exception\n");
    // test exception
    printf("test exception for param 1, uplo\n");
    SYSV_RK(&wronguplo, &N, &NRHS, A, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 2, N\n");
    SYSV_RK(&uplo, &neg, &NRHS, A, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 3, NRHS\n");
    SYSV_RK(&uplo, &N, &neg, A, &N, E, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 5, lda\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &neg, E, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 9, ldb\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &neg, W, &N, &info);
    printf("test exception for param 11, ldwork\n");
    SYSV_RK(&uplo, &N, &NRHS, A, &N, E, ipiv, B, &N, W, &neg, &info);

    free(A);
    free(B);
    free(E);
    free(W);
    free(ipiv);

    printf("test exception end\n");
}
