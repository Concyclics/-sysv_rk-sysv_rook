#include "creat_matrix.h"

void _CreatMatrix(char uplo,
                  int N,
                  int NRHS,
                  funcTdataType* A,
                  funcTdataType* B) {
    int row, col, i;
    if (uplo == 'U' || uplo == 'u') {
        for (row = 0; row < N; ++row) {
            for (col = 0; col < N; ++col) {
                if (col <= row) {
                    A[row * N + col] = ((rand() % 20000) / 10000 - 1)
#ifdef COMPLEX
                                       + ((rand() % 20000) / 10000 - 1) * I
#endif
#ifdef COMPLEX16
                                       + ((rand() % 20000) / 10000 - 1) * I
#endif
                        ;
                } else {
                    A[row * N + col] = 0;
                }
            }
        }
    } else {
        for (row = 0; row < N; ++row) {
            for (col = 0; col < N; ++col) {
                if (col >= row) {
                    A[row * N + col] = ((rand() % 20000) / 10000 - 1)
#ifdef COMPLEX
                                       + ((rand() % 20000) / 10000 - 1) * I
#endif
#ifdef COMPLEX16
                                       + ((rand() % 20000) / 10000 - 1) * I
#endif
                        ;
                } else {
                    A[row * N + col] = 0;
                }
            }
        }
    }
    // Creat matrix B.
    for (i = 0; i < NRHS * N; ++i) {
        B[i] = ((rand() % 20000) / 10000 - 1)
#ifdef COMPLEX
               + ((rand() % 20000) / 10000 - 1) * I
#endif
#ifdef COMPLEX16
               + ((rand() % 20000) / 10000 - 1) * I
#endif
            ;
    }
}
