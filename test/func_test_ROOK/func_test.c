/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The declaration of _FuncTest
 * Author: linhouzhong
 * Create: 2022-06-26
 */

#include "func_test.h"

double _TestWithSeed(int rand_seed,
                     int N,
                     int NRHS,
                     char uplo,
                     funcTdataType* A,
                     funcTdataType* A1,
                     funcTdataType* B1,
                     funcTdataType* B3,
                     funcTdataType* W,
                     funcTdataAccu* W1,
                     int* ipiv) 
{
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
    SYSV_ROOK(&uplo, &N, &NRHS, A1, &N, ipiv, B1, &N, W, &lwork,
              &info);  // lapack code
    _pot02_(&uplo, &N, &NRHS, A, &N, B1, &N, B3, &N, W1,
            &RESID);  // lapack code
    return RESID;
}

void _FuncTest(int times, char scale) {
    srand((unsigned)time(NULL));

    int size = sizeof(funcTdataType);

    for (int i = 0; i < times; ++i) {
        int N, NRHS;
        funcTdataType *A, *A1, *A2, *A3, *B1, *B2, *B3, *B4, *W;
        funcTdataAccu* W1;
        int *ipiv1, *ipiv2, info1, info2;
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

        NRHS = rand() % N + 1;  // NRHS=[1,N]

        if (((A = (funcTdataType*)malloc(size * N * N)) == NULL) ||
            ((A1 = (funcTdataType*)malloc(size * N * N)) == NULL) ||
            ((A2 = (funcTdataType*)malloc(size * N * N)) == NULL) ||
            ((A3 = (funcTdataType*)malloc(size * N * N)) == NULL) ||
            ((B1 = (funcTdataType*)malloc(size * N * NRHS)) == NULL) ||
            ((B2 = (funcTdataType*)malloc(size * N * NRHS)) == NULL) ||
            ((B4 = (funcTdataType*)malloc(size * N * NRHS)) == NULL) ||
            ((B3 = (funcTdataType*)malloc(size * N * NRHS)) == NULL) ||
            ((W = (funcTdataType*)malloc(size * N * N)) == NULL) ||
            ((W1 = (funcTdataAccu*)malloc(sizeof(funcTdataAccu) * N * N)) ==
             NULL) ||
            ((ipiv1 = (int*)malloc(sizeof(int) * N)) == NULL) ||
            ((ipiv2 = (int*)malloc(sizeof(int) * N)) == NULL)) {
            printf("malloc error!\n");
            return;
        }

        printf("  UPLO = %c\n", uplo);
        _CreatMatrix(uplo, N, NRHS, A, B1);
        // copy A to A1 and A2.
        for (int j = 0; j < N * N; ++j) {
            funcTdataType temp = A[j];
            A1[j] = temp;
            A2[j] = temp;
            A3[j] = temp;
        }
        // copy B to B1 and B2.
        for (int j = 0; j < N * NRHS; ++j) {
            funcTdataType temp = B1[j];
            B2[j] = temp;
            B3[j] = temp;
            B4[j] = temp;
        }

        // calculate result by lapack code and our code
        _sysv_rook_(&uplo, &N, &NRHS, A1, &N, ipiv1, B1, &N, W, &N,
                    &info1);  // lapack code
        _pot02_(&uplo, &N, &NRHS, A, &N, B1, &N, B3, &N, W1,
                &RESID1);  // lapack code
        SYSV_ROOK(&uplo, &N, &NRHS, A2, &N, ipiv2, B2, &N, W, &N,
                  &info2);  // our code
        _pot02_(&uplo, &N, &NRHS, A3, &N, B2, &N, B4, &N, W1,
                &RESID2);  // our code

        printf("  resid1 = %f,resid2 = %f\n", RESID1, RESID2);

        if (RESID1 >= RESID2) {
            printf("  N = %d, NRHS = %d,  succeed!\n\n", N, NRHS);
        } else {
            printf("  N = %d, NRHS = %d,  NO succeed...\n\n", N, NRHS);
        }

        free(A);
        free(A1);
        free(A2);
        free(A3);
        free(B1);
        free(B2);
        free(B3);
        free(B4);
        free(W);
        free(W1);
        free(ipiv1);
        free(ipiv2);
    }
}

void nullptr_test(int N, int M) {
    srand((unsigned)time(NULL));
    char uplo = ((rand() % 2) == 1) ? 'L' : 'U';
    dataType *A, *B, *W;
    int *ipiv, info;
    struct timespec start, end;
    double difTime, result;
    float opTRF, opTRS;
    int lin0 = 0;
    if (((A = (dataType*)malloc(sizeof(dataType) * N * N)) == NULL) ||
        ((B = (dataType*)malloc(sizeof(dataType) * N * M)) == NULL) ||
        ((W = (dataType*)malloc(sizeof(dataType) * N * N)) == NULL) ||
        ((ipiv = (int*)malloc(sizeof(int) * N)) == NULL)) {
        printf("malloc error!\n");
        return;
    }
    _CreatMatrix(uplo, N, M, A, B);
    // SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, &N, W, &N, &info);	// our
    // code
    int NRHS = N;
    printf("test nullptr\n");
    // test nullptr
    printf("test nullptr for param 1, uplo\n");
    SYSV_ROOK(NULL, &N, &NRHS, A, &N, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 2, N\n");
    SYSV_ROOK(&uplo, NULL, &NRHS, A, &N, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 3, NRHS\n");
    SYSV_ROOK(&uplo, &N, NULL, A, &N, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 4, A\n");
    SYSV_ROOK(&uplo, &N, &NRHS, NULL, &N, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 5, lda\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, NULL, ipiv, B, &N, W, &N, &info);
    printf("test nullptr for param 6, ipiv\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, NULL, B, &N, W, &N, &info);
    printf("test nullptr for param 7, B\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, NULL, &N, W, &N, &info);
    printf("test nullptr for param 8, ldb\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, NULL, W, &N, &info);
    printf("test nullptr for param 9, W\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, &N, NULL, &N, &info);
    printf("test nullptr for param 10, ldwork\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, &N, W, NULL, &info);
    printf("test nullptr for param 11, info\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, &N, W, &N, NULL);

    free(A);
    free(B);
    free(W);
    free(ipiv);
    printf("test nullptr end\n");
}

void exception_test(int N, int M) {
    srand((unsigned)time(NULL));
    char uplo = ((rand() % 2) == 1) ? 'L' : 'U';
    funcTdataType *A, *B, *W;
    int *ipiv, info;
    if (((A = (funcTdataType*)malloc(sizeof(funcTdataType) * N * N)) == NULL) ||
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
    int neg = -1;
    char wronguplo = 'a';
    printf("test exception\n");
    // test exception
    printf("test exception for param 1, uplo\n");
    SYSV_ROOK(&wronguplo, &N, &NRHS, A, &N, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 2, N\n");
    SYSV_ROOK(&uplo, &neg, &NRHS, A, &N, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 3, NRHS\n");
    SYSV_ROOK(&uplo, &N, &neg, A, &N, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 5, lda\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &neg, ipiv, B, &N, W, &N, &info);
    printf("test exception for param 8, ldb\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, &neg, W, &N, &info);
    printf("test exception for param 10, ldwork\n");
    SYSV_ROOK(&uplo, &N, &NRHS, A, &N, ipiv, B, &N, W, &neg, &info);
    
    free(A);
    free(B);
    free(W);
    free(ipiv);

    printf("test exception end\n");
}