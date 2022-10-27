#!/bin/bash
set -e

#bash build.sh 1
#bash build_test.sh


#export ASAN_OPTIONS=halt_on_error=0;log_path=/tmp/asan.log;detect_leaks=1

#export LD_PRELOAD=/usr/lib/gcc/aarch68-linux-gnu/10.3.1/libasan.so

export OMP_NUM_THREADS=32
#./test/func_test_RK/generate_func_test_RK_SINGLE.o S 10 10 U
#./test/func_test_RK/func_test_RK_SINGLE.o S 10 10 U
#./test/func_test_RK/generate_func_test_RK_DOUBLE.o S 10 10 U
#./test/func_test_RK/func_test_RK_DOUBLE.o S 10 10 U
:<<BLOCK
echo RK U
for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
./test/func_test_RK/generate_func_test_RK_$type.o S 19 10 U
./test/func_test_RK/func_test_RK_$type.o S 19 10 U
done
echo ROOK U
for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
./test/func_test_ROOK/generate_func_test_ROOK_$type.o S 19 10 U
./test/func_test_ROOK/func_test_ROOK_$type.o S 19 10 U
done

echo RK L
for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
./test/func_test_RK/generate_func_test_RK_$type.o S 19 10 L
./test/func_test_RK/func_test_RK_$type.o S 19 10 L
done
echo ROOK L
for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
    do
    ./test/func_test_ROOK/generate_func_test_ROOK_$type.o S 4 4 L
    ./test/func_test_ROOK/func_test_ROOK_$type.o S 4 4 L
    done

BLOCK

for i in {1..30..1}
do
    for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
    do
    ./test/func_test_ROOK/generate_func_test_ROOK_$type.o S 19 10 U
    ./test/func_test_ROOK/func_test_ROOK_$type.o S 19 10 U
    done
    for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
    do
    ./test/func_test_ROOK/generate_func_test_ROOK_$type.o S 19 10 L
    ./test/func_test_ROOK/func_test_ROOK_$type.o S 19 10 L
    done
done

#./test/perf_test_RK/perf_test_RK_SINGLE.o 9000 9000 U
#./test/perf_test_ROOK/perf_test_ROOK_SINGLE.o 300 300 U
#./test/perf_test_RK/perf_test_RK_SINGLE.o 5000 5000 L
#./test/perf_test_RK/perf_test_RK_SINGLE.o  L
#./run_func_test.sh
#./run_perf_test.sh

