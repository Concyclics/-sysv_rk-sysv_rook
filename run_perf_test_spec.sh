type=$1
func=$2
N=$3
NRHS=$4
UPLO=$5
numThread=$6
#tip:N>=NRHS
export OMP_NUM_THREADS=$numThread
echo NUM_THREAD $numThread
./test/perf_test_${func}/perf_test_${func}_${type}.o $N $NRHS $UPLO

