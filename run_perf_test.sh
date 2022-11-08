set -e

bash build.sh 1
bash build_test.sh

#export ASAN_OPTIONS=halt_on_error=0;log_path=/tmp/asan.log;detect_leaks=1

#export LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/10.3.1/libasan.so
echo "Running performance tests"
    echo "test RK"
for numThread in {1,48,96}
do
    export OMP_NUM_THREADS=$numThread
    echo NUM_THREAD $numThread
    echo 'test SINGLE'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_RK/perf_test_RK_SINGLE.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_RK/perf_test_RK_SINGLE.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_RK/perf_test_RK_SINGLE.o $i $i U
    done

    echo 'test DOUBLE'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_RK/perf_test_RK_DOUBLE.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_RK/perf_test_RK_DOUBLE.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_RK/perf_test_RK_DOUBLE.o $i $i U
    done

    echo 'test COMPLEX'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_RK/perf_test_RK_COMPLEX.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_RK/perf_test_RK_COMPLEX.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_RK/perf_test_RK_COMPLEX.o $i $i U
    done

    echo 'test COMPLEX16'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_RK/perf_test_RK_COMPLEX16.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_RK/perf_test_RK_COMPLEX16.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_RK/perf_test_RK_COMPLEX16.o $i $i U
    done
done


echo 'test ROOK'

for numThread in {1,48,96}
do
    export OMP_NUM_THREADS=$numThread
    echo NUM_THREAD $numThread
    echo 'test SINGLE'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_ROOK/perf_test_ROOK_SINGLE.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_ROOK/perf_test_ROOK_SINGLE.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_ROOK/perf_test_ROOK_SINGLE.o $i $i U
    done

    echo 'test DOUBLE'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_ROOK/perf_test_ROOK_DOUBLE.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_ROOK/perf_test_ROOK_DOUBLE.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_ROOK/perf_test_ROOK_DOUBLE.o $i $i U
    done

    echo 'test COMPLEX'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_ROOK/perf_test_ROOK_COMPLEX.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_ROOK/perf_test_ROOK_COMPLEX.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_ROOK/perf_test_ROOK_COMPLEX.o $i $i U
    done

    echo 'test COMPLEX16'
    echo 'test small'
    for i in $(seq 100 100 900); do
    ./test/perf_test_ROOK/perf_test_ROOK_COMPLEX16.o $i $i U
    done
    echo 'test middle'
    for i in $(seq 1000 1000 9000); do
    ./test/perf_test_ROOK/perf_test_ROOK_COMPLEX16.o $i $i U
    done
    echo 'test big'
    for i in $(seq 10000 10000 40000); do
    ./test/perf_test_ROOK/perf_test_ROOK_COMPLEX16.o $i $i U
    done
done