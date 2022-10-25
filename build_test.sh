#!/bin/bash
set -e

#export LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/10.3.1/libasan.so

cd test/perf_test_ROOK
gfortran *.c *.f -o perf_test_ROOK_SINGLE.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D SINGLE -ftrapv -fstack-protector-all
gfortran *.c *.f -o perf_test_ROOK_DOUBLE.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D DOUBLE -ftrapv -fstack-protector-all
gfortran *.c *.f -o perf_test_ROOK_COMPLEX.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D COMPLEX -ftrapv -fstack-protector-all
gfortran *.c *.f -o perf_test_ROOK_COMPLEX16.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D COMPLEX16 -ftrapv -fstack-protector-all

cd ../perf_test_RK
gfortran *.c *.f -o perf_test_RK_SINGLE.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D SINGLE -ftrapv -fstack-protector-all
gfortran *.c *.f -o perf_test_RK_DOUBLE.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D DOUBLE -ftrapv -fstack-protector-all
gfortran *.c *.f -o perf_test_RK_COMPLEX.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D COMPLEX -ftrapv -fstack-protector-all
gfortran *.c *.f -o perf_test_RK_COMPLEX16.o -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -Ofast -D COMPLEX16 -ftrapv -fstack-protector-all

cd ../func_test_RK
gfortran *.c *.f -o generate_func_test_RK_SINGLE.o -fopenmp -llapack -lblas -lkblas -lkservice -Ofast -D SINGLE -D GENERATE_DATA
gfortran *.c *.f -o func_test_RK_SINGLE.o -fopenmp -L ../.. -lSYSV -llapack -lblas -lkblas -lkservice -Ofast -D SINGLE -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs
gfortran *.c *.f -o generate_func_test_RK_DOUBLE.o -fopenmp -L ../.. -llapack -lblas -lkblas -lkservice -Ofast -D DOUBLE -D GENERATE_DATA
gfortran *.c *.f -o func_test_RK_DOUBLE.o -fopenmp -L ../.. -lSYSV -llapack -lblas -lkblas -lkservice -Ofast -D DOUBLE -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs
gfortran *.c *.f -o generate_func_test_RK_COMPLEX.o -fopenmp -L ../.. -llapack -lblas -lkblas -lkservice -Ofast -D COMPLEX -D GENERATE_DATA
gfortran *.c *.f -o func_test_RK_COMPLEX.o -fopenmp -L ../.. -lSYSV -llapack -lblas -lkblas -lkservice -Ofast -D COMPLEX -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs
gfortran *.c *.f -o generate_func_test_RK_COMPLEX16.o -fopenmp -L ../.. -llapack -lblas -lkblas -lkservice -Ofast -D COMPLEX16 -D GENERATE_DATA
gfortran *.c *.f -o func_test_RK_COMPLEX16.o -fopenmp -L ../.. -lSYSV -llapack -lblas -lkblas -lkservice -Ofast -D COMPLEX16 -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs

cd ../func_test_ROOK
gfortran *.c *.f -o generate_func_test_ROOK_SINGLE.o -fopenmp -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D SINGLE -D GENERATE_DATA
gfortran *.c *.f -o func_test_ROOK_SINGLE.o -fopenmp -L ../.. -lSYSV -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D SINGLE -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs
gfortran *.c *.f -o generate_func_test_ROOK_DOUBLE.o -fopenmp -L ../.. -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D DOUBLE -D GENERATE_DATA
gfortran *.c *.f -o func_test_ROOK_DOUBLE.o -fopenmp -L ../.. -lSYSV -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D DOUBLE -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs
gfortran *.c *.f -o generate_func_test_ROOK_COMPLEX.o -fopenmp -L ../.. -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D COMPLEX -D GENERATE_DATA
gfortran *.c *.f -o func_test_ROOK_COMPLEX.o -fopenmp -L ../.. -lSYSV -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D COMPLEX -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs
gfortran *.c *.f -o generate_func_test_ROOK_COMPLEX16.o -fopenmp -L ../.. -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D COMPLEX16 -D GENERATE_DATA
gfortran *.c *.f -o func_test_ROOK_COMPLEX16.o -fopenmp -L ../.. -lSYSV -llapack -lcblas -lblas -lkblas -lkservice -Ofast -D COMPLEX16 -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs