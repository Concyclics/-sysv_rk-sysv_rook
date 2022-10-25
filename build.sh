#!/bin/bash
set -e

if (($1 == 1))
then 

flags="-fPIC -O2 -fopenmp -ftrapv -fstack-protector-all -I ../../include/SYSV_RK"

cd src/SYSV_RK
#生成.o文件
gcc -c ilaenv.c -o ilaenv.o $flags  
gcc -c LASYF_RK.c -o LASYF_RK_SINGLE.o $flags -D SINGLE  
gcc -c SYSV_RK.c -o SYSV_RK_SINGLE.o $flags -D SINGLE  
gcc -c SYTF2_RK.c -o SYTF2_RK_SINGLE.o $flags -D SINGLE  
gcc -c SYTRF_RK.c -o SYTRF_RK_SINGLE.o $flags -D SINGLE  
gcc -c SYTRS_3.c -o SYTRS_3_SINGLE.o $flags -D SINGLE  

gcc -c LASYF_RK.c -o LASYF_RK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYSV_RK.c -o SYSV_RK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYTF2_RK.c -o SYTF2_RK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYTRF_RK.c -o SYTRF_RK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYTRS_3.c -o SYTRS_3_DOUBLE.o $flags -D DOUBLE  

gcc -c LASYF_RK.c -o LASYF_RK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYSV_RK.c -o SYSV_RK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYTF2_RK.c -o SYTF2_RK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYTRF_RK.c -o SYTRF_RK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYTRS_3.c -o SYTRS_3_COMPLEX.o $flags -D COMPLEX  

gcc -c LASYF_RK.c -o LASYF_RK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYSV_RK.c -o SYSV_RK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYTF2_RK.c -o SYTF2_RK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYTRF_RK.c -o SYTRF_RK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYTRS_3.c -o SYTRS_3_COMPLEX16.o $flags -D COMPLEX16  


flags="-fPIC -O2 -fopenmp -ftrapv -fstack-protector-all -I ../../include/SYSV_ROOK"
cd ../SYSV_ROOK
#生成.o文件
gcc -c LASYF_ROOK.c -o LASYF_ROOK_SINGLE.o $flags -D SINGLE  
gcc -c SYSV_ROOK.c -o SYSV_ROOK_SINGLE.o $flags -D SINGLE  
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_SINGLE.o $flags -D SINGLE  
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_SINGLE.o $flags -D SINGLE  
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_SINGLE.o $flags -D SINGLE  

gcc -c LASYF_ROOK.c -o LASYF_ROOK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYSV_ROOK.c -o SYSV_ROOK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_DOUBLE.o $flags -D DOUBLE  
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_DOUBLE.o $flags -D DOUBLE  

gcc -c LASYF_ROOK.c -o LASYF_ROOK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYSV_ROOK.c -o SYSV_ROOK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_COMPLEX.o $flags -D COMPLEX  
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_COMPLEX.o $flags -D COMPLEX  

gcc -c LASYF_ROOK.c -o LASYF_ROOK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYSV_ROOK.c -o SYSV_ROOK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_COMPLEX16.o $flags -D COMPLEX16  
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_COMPLEX16.o $flags -D COMPLEX16  

#生成静态库、动态库
ar -rsv ../../libSYSV.a *.o ../SYSV_RK/*.o
gcc -shared -o ../../libSYSV.so *.o ../SYSV_RK/*.o $flags

else

flags="-fPIC  -O2 -fopenmp -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs -I ../../include/SYSV_RK"

cd src/SYSV_RK
#生成.o文件
gcc -c ilaenv.c -o ilaenv.o $flags 
gcc -c LASYF_RK.c -o LASYF_RK_SINGLE.o $flags -D SINGLE
gcc -c SYSV_RK.c -o SYSV_RK_SINGLE.o $flags -D SINGLE
gcc -c SYTF2_RK.c -o SYTF2_RK_SINGLE.o $flags -D SINGLE
gcc -c SYTRF_RK.c -o SYTRF_RK_SINGLE.o $flags -D SINGLE
gcc -c SYTRS_3.c -o SYTRS_3_SINGLE.o $flags -D SINGLE

gcc -c LASYF_RK.c -o LASYF_RK_DOUBLE.o $flags -D DOUBLE
gcc -c SYSV_RK.c -o SYSV_RK_DOUBLE.o $flags -D DOUBLE
gcc -c SYTF2_RK.c -o SYTF2_RK_DOUBLE.o $flags -D DOUBLE
gcc -c SYTRF_RK.c -o SYTRF_RK_DOUBLE.o $flags -D DOUBLE
gcc -c SYTRS_3.c -o SYTRS_3_DOUBLE.o $flags -D DOUBLE

gcc -c LASYF_RK.c -o LASYF_RK_COMPLEX.o $flags -D COMPLEX
gcc -c SYSV_RK.c -o SYSV_RK_COMPLEX.o $flags -D COMPLEX
gcc -c SYTF2_RK.c -o SYTF2_RK_COMPLEX.o $flags -D COMPLEX
gcc -c SYTRF_RK.c -o SYTRF_RK_COMPLEX.o $flags -D COMPLEX
gcc -c SYTRS_3.c -o SYTRS_3_COMPLEX.o $flags -D COMPLEX

gcc -c LASYF_RK.c -o LASYF_RK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYSV_RK.c -o SYSV_RK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYTF2_RK.c -o SYTF2_RK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYTRF_RK.c -o SYTRF_RK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYTRS_3.c -o SYTRS_3_COMPLEX16.o $flags -D COMPLEX16


flags="-fPIC  -O2 -fopenmp -ftrapv -fstack-protector-all -fno-omit-frame-pointer -fsanitize=address -ftest-coverage -fprofile-arcs -I ../../include/SYSV_ROOK"
cd ../SYSV_ROOK
#生成.o文件
gcc -c LASYF_ROOK.c -o LASYF_ROOK_SINGLE.o $flags -D SINGLE
gcc -c SYSV_ROOK.c -o SYSV_ROOK_SINGLE.o $flags -D SINGLE
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_SINGLE.o $flags -D SINGLE
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_SINGLE.o $flags -D SINGLE
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_SINGLE.o $flags -D SINGLE

gcc -c LASYF_ROOK.c -o LASYF_ROOK_DOUBLE.o $flags -D DOUBLE
gcc -c SYSV_ROOK.c -o SYSV_ROOK_DOUBLE.o $flags -D DOUBLE
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_DOUBLE.o $flags -D DOUBLE
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_DOUBLE.o $flags -D DOUBLE
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_DOUBLE.o $flags -D DOUBLE

gcc -c LASYF_ROOK.c -o LASYF_ROOK_COMPLEX.o $flags -D COMPLEX
gcc -c SYSV_ROOK.c -o SYSV_ROOK_COMPLEX.o $flags -D COMPLEX
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_COMPLEX.o $flags -D COMPLEX
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_COMPLEX.o $flags -D COMPLEX
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_COMPLEX.o $flags -D COMPLEX

gcc -c LASYF_ROOK.c -o LASYF_ROOK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYSV_ROOK.c -o SYSV_ROOK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYTF2_ROOK.c -o SYTF2_ROOK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYTRF_ROOK.c -o SYTRF_ROOK_COMPLEX16.o $flags -D COMPLEX16
gcc -c SYTRS_ROOK.c -o SYTRS_ROOK_COMPLEX16.o $flags -D COMPLEX16

#生成静态库、动态库
ar -rsv ../../libSYSV.a *.o ../SYSV_RK/*.o
gcc -shared -fPIC -o ../../libSYSV.so *.o ../SYSV_RK/*.o $flags

fi

