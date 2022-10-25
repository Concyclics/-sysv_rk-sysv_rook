/*
gfortran -c _opla.f -o _opla.o
gfortran main.c perf_test.c creat_matrix.c _opla.o SYSV_RK.c SYTF2_RK.c
SYTRF_RK.c SYTRS_3.c LASYF_RK.c ilaenv.c -o main -L /usr/local/kml/lib
-lkservice -L /usr/local/lapacke/lib -llapack -L /usr/local/kml/lib/kblas/omp
-lkblas -fopenmp -D
*/
#include "perf_test.h"
int main(int argc, char* argv[]) {
    int N, NRHS;
    char UPLO;
    srand((unsigned)time(NULL));
    N = atoi(argv[1]);
    NRHS = atoi(argv[2]);
    UPLO = argv[3][0];
    printf("\nN=%d,NRHS=%d,UPLO=%c\n", N, NRHS, UPLO);
    perfTest(N, NRHS, UPLO);
    return 0;
}