#include "perf_test.h"
//#include <gperftools/profiler.h>
void perfTest(int N, int NRHS, char UPLO) {
    // srand((unsigned)time(NULL));
    // int n, nrhs;
    double gflops;
    gflops = calGflops(N, NRHS, UPLO);
    printf("gflops=%lf\n", gflops);
}

double calGflops(int N, int M, char uplo) {
    dataType *A, *B, *W;
    int *ipiv, info;
    char TRF[] = SYTRF_NAME, TRS[] = SYTRS_NAME;
    struct timeval tvStart, tvEnd;
    struct timezone tzStart, tzEnd;
    double result;
    double opTRF, opTRS;
    int lin0 = 0;
    if (((A = (dataType*)malloc(sizeof(dataType) * N * N)) == NULL) ||
        ((B = (dataType*)malloc(sizeof(dataType) * N * M)) == NULL) ||
        ((W = (dataType*)malloc(sizeof(dataType) * 64 * N)) == NULL) ||
        ((ipiv = (int*)malloc(sizeof(int) * N)) == NULL)) {
        printf("malloc error!\n");
        return -1;
    }
    _CreatMatrix(uplo, N, M, A, B);
    int lwork = N * 64;
    // ProfilerStart("perf_test.prof");
    gettimeofday(&tvStart, &tzStart);
    SYSV_ROOK(&uplo, &N, &M, A, &N, ipiv, B, &N, W, &lwork, &info);
    gettimeofday(&tvEnd, &tzEnd);
    // ProfilerStop();
    opTRF = _opla_(TRF, &N, &M, &lin0, &lin0, &lin0);
    opTRS = _opla_(TRS, &N, &M, &lin0, &lin0, &lin0);
    result = (double)(opTRF + opTRS) /
             ((double)(tvEnd.tv_sec - tvStart.tv_sec) * (double)1000000000 +
              (double)(tvEnd.tv_usec - tvStart.tv_usec) * (double)1000);

    free(A);
    free(B);
    free(W);
    free(ipiv);
    return result;
}
