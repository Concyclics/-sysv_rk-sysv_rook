#!/bin/bash
set -e

if (($1 == 1))
then 
selectFlag="-ftrapv -fstack-protector-all"
else
selectFlag="-ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs"
fi

#静态库
flags="-O2 -fopenmp -ftrapv $selectFlag -I ../../include/SYSV_RK -I ../../include/kml"
cd src/SYSV_RK
gcc -c ilaenv.c -o ilaenv.o $flags  
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
gcc -c LASYF_RK.c -o LASYF_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYSV_RK.c -o SYSV_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYTF2_RK.c -o SYTF2_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRF_RK.c -o SYTRF_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRS_3.c -o SYTRS_3_$DFlag.o $flags -D $DFlag  
done

flags="-O2 -fopenmp -ftrapv $selectFlag -I ../../include/SYSV_ROOK -I ../../include/kml"
cd ../SYSV_ROOK
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
gcc -c LASYF_ROOK.c -o LASYF_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYSV_ROOK.c -o SYSV_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_$DFlag.o $flags -D $DFlag  
done 

#生成静态库
ar -rsv ../../libSYSV.a *.o ../SYSV_RK/*.o

#动态库
flags="-fPIC -O2 -fopenmp -ftrapv $selectFlag -I ../../include/SYSV_RK -I ../../include/kml"
cd ../SYSV_RK
gcc -c ilaenv.c -o ilaenv.o $flags  
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
gcc -c LASYF_RK.c -o LASYF_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYSV_RK.c -o SYSV_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYTF2_RK.c -o SYTF2_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRF_RK.c -o SYTRF_RK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRS_3.c -o SYTRS_3_$DFlag.o $flags -D $DFlag  
done


flags="-fPIC -O2 -fopenmp -ftrapv $selectFlag -I ../../include/SYSV_ROOK -I ../../include/kml"
cd ../SYSV_ROOK
for DFlag in {'SINGLE','DOUBLE','COMPLEX','COMPLEX16'}
do
gcc -c LASYF_ROOK.c -o LASYF_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYSV_ROOK.c -o SYSV_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_$DFlag.o $flags -D $DFlag  
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_$DFlag.o $flags -D $DFlag  
done 

#生成动态库
gcc -shared -o ../../libSYSV.so *.o ../SYSV_RK/*.o $flags


