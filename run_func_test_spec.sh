type=$1
func=$2
N=$3
NRHS=$4
UPLO=$5
#tip:N>=NRHS
./test/func_test_${func}/generate_func_test_${func}_${type}.o S $N $NRHS $UPLO
./test/func_test_${func}/func_test_${func}_${type}.o S $N $NRHS $UPLO

