set -e

export ASAN_OPTIONS=halt_on_error=0;log_path=/tmp/asan.log;detect_leaks=1

export LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/10.3.1/libasan.so

echo "Running functional tests"

echo "test RK"

echo 'test SINGLE'
./test/func_test_RK/func_test_RK_SINGLE.o N
./test/func_test_RK/func_test_RK_SINGLE.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_RK/generate_func_test_RK_SINGLE.o S $i $i U
./test/func_test_RK/func_test_RK_SINGLE.o S $i $i U
./test/func_test_RK/generate_func_test_RK_SINGLE.o S $i $i L
./test/func_test_RK/func_test_RK_SINGLE.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_RK/generate_func_test_RK_SINGLE.o S $i $i U
./test/func_test_RK/func_test_RK_SINGLE.o S $i $i U
./test/func_test_RK/generate_func_test_RK_SINGLE.o S $i $i L
./test/func_test_RK/func_test_RK_SINGLE.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_RK/generate_func_test_RK_SINGLE.o S $i $i U
./test/func_test_RK/func_test_RK_SINGLE.o S $i $i U
./test/func_test_RK/generate_func_test_RK_SINGLE.o S $i $i L
./test/func_test_RK/func_test_RK_SINGLE.o S $i $i L
done

echo 'test DOUBLE'
./test/func_test_RK/func_test_RK_DOUBLE.o N
./test/func_test_RK/func_test_RK_DOUBLE.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S $i $i U
./test/func_test_RK/func_test_RK_DOUBLE.o S $i $i U
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S $i $i L
./test/func_test_RK/func_test_RK_DOUBLE.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S $i $i U
./test/func_test_RK/func_test_RK_DOUBLE.o S $i $i U
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S $i $i L
./test/func_test_RK/func_test_RK_DOUBLE.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S $i $i U
./test/func_test_RK/func_test_RK_DOUBLE.o S $i $i U
./test/func_test_RK/generate_func_test_RK_DOUBLE.o S $i $i L
./test/func_test_RK/func_test_RK_DOUBLE.o S $i $i L
done

echo 'test COMPLEX'
./test/func_test_RK/func_test_RK_COMPLEX.o N
./test/func_test_RK/func_test_RK_COMPLEX.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_RK/generate_func_test_RK_COMPLEX.o S $i $i U
./test/func_test_RK/func_test_RK_COMPLEX.o S $i $i U
./test/func_test_RK/generate_func_test_RK_COMPLEX.o S $i $i L
./test/func_test_RK/func_test_RK_COMPLEX.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_RK/generate_func_test_RK_COMPLEX.o S $i $i U
./test/func_test_RK/func_test_RK_COMPLEX.o S $i $i U
./test/func_test_RK/generate_func_test_RK_COMPLEX.o S $i $i L
./test/func_test_RK/func_test_RK_COMPLEX.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_RK/generate_func_test_RK_COMPLEX.o S $i $i U
./test/func_test_RK/func_test_RK_COMPLEX.o S $i $i U
./test/func_test_RK/generate_func_test_RK_COMPLEX.o S $i $i L
./test/func_test_RK/func_test_RK_COMPLEX.o S $i $i L
done

echo 'test COMPLEX16'
./test/func_test_RK/func_test_RK_COMPLEX16.o N
./test/func_test_RK/func_test_RK_COMPLEX16.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_RK/generate_func_test_RK_COMPLEX16.o S $i $i U
./test/func_test_RK/func_test_RK_COMPLEX16.o S $i $i U
./test/func_test_RK/generate_func_test_RK_COMPLEX16.o S $i $i L
./test/func_test_RK/func_test_RK_COMPLEX16.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_RK/generate_func_test_RK_COMPLEX16.o S $i $i U
./test/func_test_RK/func_test_RK_COMPLEX16.o S $i $i U
./test/func_test_RK/generate_func_test_RK_COMPLEX16.o S $i $i L
./test/func_test_RK/func_test_RK_COMPLEX16.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_RK/generate_func_test_RK_COMPLEX16.o S $i $i U
./test/func_test_RK/func_test_RK_COMPLEX16.o S $i $i U
./test/func_test_RK/generate_func_test_RK_COMPLEX16.o S $i $i L
./test/func_test_RK/func_test_RK_COMPLEX16.o S $i $i L
done

echo 'test ROOK'

echo 'test SINGLE'
./test/func_test_ROOK/func_test_ROOK_SINGLE.o N
./test/func_test_ROOK/func_test_ROOK_SINGLE.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_ROOK/generate_func_test_ROOK_SINGLE.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_SINGLE.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_SINGLE.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_SINGLE.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_ROOK/generate_func_test_ROOK_SINGLE.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_SINGLE.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_SINGLE.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_SINGLE.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_ROOK/generate_func_test_ROOK_SINGLE.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_SINGLE.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_SINGLE.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_SINGLE.o S $i $i L
done

echo 'test DOUBLE'
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o N
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_ROOK/generate_func_test_ROOK_DOUBLE.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_DOUBLE.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_ROOK/generate_func_test_ROOK_DOUBLE.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_DOUBLE.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_ROOK/generate_func_test_ROOK_DOUBLE.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_DOUBLE.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_DOUBLE.o S $i $i L
done

echo 'test COMPLEX'
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o N
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_COMPLEX.o S $i $i L
done

echo 'test COMPLEX16'
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o N
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o E
echo 'test small'
for i in $(seq 100 100 900); do
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX16.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX16.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o S $i $i L
done
echo 'test medium'
for i in $(seq 1000 1000 9000); do
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX16.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX16.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o S $i $i L
done
echo 'test large'
for i in $(seq 10000 10000 40000); do
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX16.o S $i $i U
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o S $i $i U
./test/func_test_ROOK/generate_func_test_ROOK_COMPLEX16.o S $i $i L
./test/func_test_ROOK/func_test_ROOK_COMPLEX16.o S $i $i L
done