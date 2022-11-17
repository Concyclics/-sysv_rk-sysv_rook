#!/bin/bash
set -e

bash build.sh 2
bash build_test.sh

export ASAN_OPTIONS=log_path=./asan.log:detect_leaks=1
export OMP_NUM_THREADS=32
#export LD_PRELOAD=/usr/lib/gcc/aarch64-linux-gnu/10.3.1/libasan.so

echo "Running functional tests"

echo "test RK"
for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
    echo "test $type"
    ./test/func_test_RK/func_test_RK_$type.o N
    ./test/func_test_RK/func_test_RK_$type.o E
    ./test/func_test_RK/func_test_RK_$type.o C
:<<BLOCK
    echo 'test small'
    ((time=$RANDOM%11+20))
    echo "test time = $time"
    for ((i=0;i<$time;i++)); do
        ((UPLO=$RANDOM%2+1))
        if [ $UPLO -eq 1 ]
        then
            UPLO=L
        else
            UPLO=U
        fi
        
        ((N=$RANDOM%900+100))
        ((NRHS=$RANDOM%($N-100)+100))
        echo N=$N
        echo NRHS=$NRHS
        echo UPLO=$UPLO
        ./test/func_test_RK/generate_func_test_RK_$type.o S $N $NRHS $UPLO
        ./test/func_test_RK/func_test_RK_$type.o S $N $NRHS $UPLO
    done

    echo 'test medium'
    ((time=$RANDOM%11+10))
    echo "test time = $time"
    for ((i=0;i<$time;i++)); do
        ((UPLO=$RANDOM%2+1))
        if [ $UPLO -eq 1 ]
        then
            UPLO=L
        else
            UPLO=U
        fi
        ((N=$RANDOM%9000+1000))
        ((NRHS=$RANDOM%($N-1000)+1000))
        echo N=$N
        echo NRHS=$NRHS
        echo UPLO=$UPLO
        ./test/func_test_RK/generate_func_test_RK_$type.o S $N $NRHS $UPLO
        ./test/func_test_RK/func_test_RK_$type.o S $N $NRHS $UPLO
    done

    echo 'test large'
    ((time=$RANDOM%5+1))
    echo "test time = $time"
    for ((i=0;i<$time;i++)); do
        ((UPLO=$RANDOM%2+1))
        if [ $UPLO -eq 1 ]
        then
            UPLO=L
        else
            UPLO=U
        fi
        ((N=$RANDOM%30000+10000))
        ((NRHS=$RANDOM%($N-10000)+10000))
        echo N=$N
        echo NRHS=$NRHS
        echo UPLO=$UPLO
        ./test/func_test_RK/generate_func_test_RK_$type.o S $N $NRHS $UPLO
        ./test/func_test_RK/func_test_RK_$type.o S $N $NRHS $UPLO
    done
BLOCK
done


echo 'test ROOK'

for type in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
    echo "test $type"
    ./test/func_test_ROOK/func_test_ROOK_$type.o N
    ./test/func_test_ROOK/func_test_ROOK_$type.o E
    ./test/func_test_ROOK/func_test_ROOK_$type.o C
:<<BLOCK
    echo 'test small'
    ((time=$RANDOM%11+20))
    echo "test time = $time"
    for ((i=0;i<$time;i++)); do
        ((UPLO=$RANDOM%2+1))
        if [ $UPLO -eq 1 ]
        then
            UPLO=L
        else
            UPLO=U
        fi
        ((N=$RANDOM%900+100))
        ((NRHS=$RANDOM%($N-100)+100))
        echo N=$N
        echo NRHS=$NRHS
        echo UPLO=$UPLO
        ./test/func_test_ROOK/generate_func_test_ROOK_$type.o S $N $NRHS $UPLO
        ./test/func_test_ROOK/func_test_ROOK_$type.o S $N $NRHS $UPLO
    done

    echo 'test medium'
    ((time=$RANDOM%11+10))
    echo "test time = $time"
    for ((i=0;i<$time;i++)); do
        ((UPLO=$RANDOM%2+1))
        if [ $UPLO -eq 1 ]
        then
            UPLO=L
        else
            UPLO=U
        fi
        ((N=$RANDOM%9000+1000))
        ((NRHS=$RANDOM%($N-1000)+1000))
        echo N=$N
        echo NRHS=$NRHS
        echo UPLO=$UPLO
        ./test/func_test_ROOK/generate_func_test_ROOK_$type.o S $N $NRHS $UPLO
        ./test/func_test_ROOK/func_test_ROOK_$type.o S $N $NRHS $UPLO
    done

    echo 'test large'
    ((time=$RANDOM%5+1))
    echo "test time = $time"
    echo "test time = $time"
    for ((i=0;i<$time;i++)); do
        ((UPLO=$RANDOM%2+1))
        if [ $UPLO -eq 1 ]
        then
            UPLO=L
        else
            UPLO=U
        fi
        ((N=$RANDOM%30000+10000))
        ((NRHS=$RANDOM%($N-10000)+10000))
        echo N=$N
        echo NRHS=$NRHS
        echo UPLO=$UPLO
        ./test/func_test_ROOK/generate_func_test_ROOK_$type.o S $N $NRHS $UPLO
        ./test/func_test_ROOK/func_test_ROOK_$type.o S $N $NRHS $UPLO
    done
BLOCK
done
