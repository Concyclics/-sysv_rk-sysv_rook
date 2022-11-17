/*
 * Copyright (c) SCUT. 2022. All rights reserved.
 * Description: The main of function test.
 * Author: linhouzhong
 * Create: 2022-06-26
 */

// gfortran main.c func_test.c creat_matrix.c _pot02.o SYSV_RK.c SYTF2_RK.c
// SYTRF_RK.c SYTRS_3.c LASYF_RK.c ilaenv.c   -o main -L /usr/local/lapacke/lib
// -llapack  -L /usr/local/cblas -lcblas -lblas -L /usr/local/kml/lib -lkservice
// -D SINGLE -fopenmp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "func_test.h"

int main(int argc, char* argv[]) {
    char method = argv[1][0];
    if (method == 'H' || method == 'h') {
        printf("func_test help:\n");
        printf("N: test nullptr\n");
        printf("E: test exception\n");
        printf("S: solve matrix\n with 3 args: N, NRHS, UPLO\n");
        return 0;
    }
    if (method == 'N' || method == 'n') {
        printf("func_test nullptr:\n");
        nullptr_test(10, 10);
        return 0;
    }
    if (method == 'E' || method == 'e') {
        printf("func_test exception:\n");
        exception_test(10, 10);
        return 0;
    }
    if (method == 'C' || method == 'c') {
        printf("cover:\n");
        return cover();
    }
    if (method != 'S' && method != 's') {
        printf("func_test error: wrong method\n");
        return 0;
    }

    int N, NRHS;
    char UPLO;
    N = atoi(argv[2]);
    NRHS = atoi(argv[3]);
    UPLO = argv[4][0];
    funcTdataType *A, *A1, *E, *B1, *B3, *W;
    funcTdataAccu* W1;
    int* ipiv;
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

#ifdef GENERATE_DATA

    int setSeed = 19942;
    // srand(time(NULL));
    srand(setSeed);
    int seed = rand();
    double result;
    freopen("seeds_answers_TestTypeName.txt", "w", stdout);
    result = _TestWithSeed(seed, N, NRHS, UPLO, A, A1, E, B1, B3, W, W1, ipiv);
    printf("%d %f\n", seed, result);
    result = _TestWithSeedGet04(seed, N, NRHS, UPLO, A, A1, E, B1, B3, W, W1,
                                ipiv);
    printf("%d %f\n", seed, result);
    fclose(stdout);
    free(A);
    free(A1);
    free(E);
    free(B1);
    free(B3);
    free(W);
    free(W1);
    free(ipiv);
#else

    freopen("seeds_answers_TestTypeName.txt", "r", stdin);
    int seed;
    double result_OS, result_US;
    scanf("%d %lf", &seed, &result_OS);
    result_US =
        _TestWithSeed(seed, N, NRHS, UPLO, A, A1, E, B1, B3, W, W1, ipiv);
    printf("test pot02\n");
    printf("test N = %d, NRHS = %d, UPLO = %c\n", N, NRHS, UPLO);
    printf("result_OS = %f, result_US = %f\n", result_OS, result_US);
    if (result_OS > result_US) {
        printf("success!\n");
    } else {
        printf("fail!\n");
    }
    scanf("%d %lf", &seed, &result_OS);
    result_US =
        _TestWithSeedGet04(seed, N, NRHS, UPLO, A, A1, E, B1, B3, W, W1, ipiv);
    printf("test get04\n");
    printf("test N = %d, NRHS = %d, UPLO = %c\n", N, NRHS, UPLO);
    printf("result_OS = %f, result_US = %f\n", result_OS, result_US);
    if (result_OS > result_US) {
        printf("success!\n");
    } else {
        printf("fail!\n");
    }


    fclose(stdin);
    free(A);
    free(A1);
    free(E);
    free(B1);
    free(B3);
    free(W);
    free(W1);
    free(ipiv);
#endif
    return 0;
}