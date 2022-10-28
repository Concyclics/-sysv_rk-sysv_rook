#!/bin/bash
set -e

#export LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/10.3.1/libasan.so
selectFlag="-ftrapv -fstack-protector-all"
flags=" -fopenmp -L ../.. -lSYSV -lkblas -llapack -lkservice -O2 $selectFlag"
cd test/perf_test_ROOK
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
gfortran *.c *.f -o perf_test_ROOK_$DFlag.o $flags -D $DFlag
done

cd ../perf_test_RK
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
gfortran *.c *.f -o perf_test_RK_$DFlag.o $flags -D $DFlag
done

selectFlag="-ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs"
flags1="-fopenmp -llapack -lblas -lkservice -O2"
flags2="-fopenmp -L ../.. -lSYSV -llapack -lkblas -lblas -lkservice -O2 $selectFlag"
cd ../func_test_RK
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do 
gfortran *.c *.f -o generate_func_test_RK_$DFlag.o $flags1 -D GENERATE_DATA -D $DFlag
gfortran *.c *.f -o func_test_RK_$DFlag.o $flags2 -D $DFlag
done

cd ../func_test_ROOK
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do 
gfortran *.c *.f -o generate_func_test_ROOK_$DFlag.o $flags1 -D GENERATE_DATA -D $DFlag
gfortran *.c *.f -o func_test_ROOK_$DFlag.o $flags2 -D $DFlag
done