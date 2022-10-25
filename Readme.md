
# SYSV使用说明
## 接口介绍
本接口由SYSV_RK以及SYSV_ROOK组成，包含float,double,float complex,double complex四种数据类型，总共8个接口。
****
## 组成
本文件主要由源文件src以及测试文件test组成。src中包含SYSV_RK以及SYSV_ROOK。test中包含四种数据类型的测试。
****
## 编译说明
1、在SYSV目录下，可以选择两种编译组合：

(1)仅加入安全选项
```bash
build.sh 1
```
即可将源文件编译成libSYSV.a以及libSYSV.so两个库。

(2)加入三个选项
```bash
build.sh 2
```
即可将源文件编译成libSYSV.a以及libSYSV.so两个库。


2、在在SYSV目录下，使用命令
```bash
build_test.sh
```
即可将编译各个测试代码。
****
## 功能测试使用说明（以RK下的SINGLE类型为例）
在SYSV目录下，使用命令
```bash
run_func_test.sh
```
其中，若想单独测试某一功能，则使用命令：（以RK下的SINGLE类型为例）
        
测试空指针
```bash
./test/func_test_RK/func_test_RK_SINGLE.o N
```
测试异常输入
```bash
./test/func_test_RK/func_test_RK_SINGLE.o E
```
测试某一规模矩阵
```bash
./test/func_test_RK/generate_func_test_RK_SINGLE.o S N NRHS UPLO
./test/func_test_RK/func_test_RK_SINGLE.o S N NRHS UPLO
```

        
****
## 性能测试使用说明
在SYSV目录下，使用命令
```bash
run_perf_test.sh
```
其中，若想单独测试某一规模，则使用命令：（以RK下的SINGLE类型为例）
```bash
./test/perf_test_RK/perf_test_RK_SINGLE.o N NRHS UPLO
```
其中N和NRHS为规模数字，UPLO为U或者L。




