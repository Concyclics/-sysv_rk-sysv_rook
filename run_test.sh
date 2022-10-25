#!/bin/bash
set -e

bash build.sh 1
bash build_test.sh


#export ASAN_OPTIONS=halt_on_error=0;log_path=/tmp/asan.log;detect_leaks=1

#export LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/10.3.1/libasan.so

export OMP_NUM_THREADS=32
#./test/func_test_RK/generate_func_test_RK_SINGLE.o S 232 232 U
#./test/func_test_RK/func_test_RK_SINGLE.o S 232 232 U
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S 232 232 U
./test/func_test_RK/func_test_RK_DOUBLE.o S 232 232 U
./test/perf_test_RK/perf_test_RK_SINGLE.o 9000 9000 U
./test/perf_test_RK/perf_test_RK_SINGLE.o 9000 9000 L
#./test/perf_test_RK/perf_test_RK_SINGLE.o 5000 5000 L
#./test/perf_test_RK/perf_test_RK_SINGLE.o  L
#./run_func_test.sh
#./run_perf_test.sh

