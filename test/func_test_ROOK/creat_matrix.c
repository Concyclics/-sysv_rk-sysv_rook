#include "creat_matrix.h"

void _CreatMatrix(char uplo,
                  int N,
                  int NRHS,
                  funcTdataType* A,
                  funcTdataType* B) {
    int row, col, i;
    // srand((unsigned)time(NULL));
    //  Creat matrix A.
    if (uplo == 'U' || uplo == 'u') {
        //#pragma omp parallel for schedule(static) collapse(2) private(row,
        //col)
        for (row = 0; row < N; ++row) {
            for (col = 0; col < N; ++col) {
                if (col <= row) {
                    A[row * N + col] =
                        ((rand() % 20000) - 10000)
#ifdef COMPLEX
                        + (((rand() % 20000) - 10000) / 100.0) * I
#endif
#ifdef COMPLEX16
                        + (((rand() % 20000) - 10000) / 100.0) * I
#endif
                        ;
                } else {
                    A[row * N + col] = 0;
                }
            }
        }
    } else {
        //#pragma omp parallel for schedule(static) collapse(2) private(row,
        //col)
        for (row = 0; row < N; ++row) {
            for (col = 0; col < N; ++col) {
                if (col >= row) {
                    A[row * N + col] =
                        ((rand() % 20000) - 10000)
#ifdef COMPLEX
                        + (((rand() % 20000) - 10000) / 100.0) * I
#endif
#ifdef COMPLEX16
                        + (((rand() % 20000) - 10000) / 100.0) * I
#endif
                        ;
                } else {
                    A[row * N + col] = 0;
                }
            }
        }
    }
    // Creat matrix B.
    //#pragma omp parallel for schedule(static) private(i)
    for (i = 0; i < NRHS * N; ++i) {
        B[i] = ((rand() % 20000) - 10000)
#ifdef COMPLEX
               + (((rand() % 20000) - 10000) / 100.0) * I
#endif
#ifdef COMPLEX16
               + (((rand() % 20000) - 10000) / 100.0) * I
#endif
            ;
    }
}
