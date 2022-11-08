
# SYSV使用说明
## 接口介绍
本接口由SYSV_RK以及SYSV_ROOK组成，包含float,double,float complex,double complex四种数据类型，总共8个接口。

| 数据类型           |      接口名称      
|----------         |:-------------:
| float             |  ssysv_rk/ssysv_rook 
| double            |  dsysv_rk/dsysv_rook  
| float complex     |  csysv_rk/csysv_rook 
| double complex    |  zsysv_rk/zsysv_rook

****
## 存储类型
普通矩阵按照Fortran的列优先方式存储（column-major），数据结构中涉及到的索引从1开始（如ipiv中保存的值）。

压缩存储矩阵（packed）：以列优先方式存储，对于实对称、共轭对称、上三角或下三角矩阵，只需一半的元素即可确定整个矩阵，因此可使用压缩格式只保存这一部分数据。对于大小的压缩存储矩阵A，对应元素的存储位置如下：

上三角存储（uplo='U'，且1<=i<=j<=n）:i-1+j*(j-1)/2

下三角存储（uplo='L'，且1<=j<=i<=n）:i-1+(j-1)*(2n-j)/2

****
## 组成
本文件主要由源文件src以及测试文件test组成。src中包含SYSV_RK以及SYSV_ROOK。test中包含四种数据类型的测试。
****
## 编译说明
1、在SYSV目录下，可以选择两种编译组合：

(1)仅加入安全选项
```
bash build.sh 1
```
即可将源文件编译成libSYSV.a以及libSYSV.so两个库。

(2)加入三个选项
```
bash build.sh 2
```
即可将源文件编译成libSYSV.a以及libSYSV.so两个库。


2、在在SYSV目录下，使用命令
```
bash build_test.sh
```
即可将编译各个测试代码。
****
## 功能测试使用说明
### 完整功能测试（均在SYSV目录下使用）
在SYSV目录下，使用以下命令，其中使用了对应的编译选项进行编译（输出中result_os为开源代码误差，result_us为开发代码误差）
```
bash run_func_test.sh
```
测试空指针
```
bash ./test/func_test_RK/func_test_RK_SINGLE.o N
```
测试异常输入
```
bash ./test/func_test_RK/func_test_RK_SINGLE.o E
```
### 单独测试某一规模
若想单独测试某一功能，则使用命令：
```
bash run_func_test_spec.sh type func N NRHS UPLO
```
以RK函数的SINGLE类型，N=100，NRHS=50，UPLO=L为例:(注意所有字母要大写)
```
bash run_func_test_spec.sh SINGLE RK 100 50 L
```
       
****
## 性能测试使用说明
### 完整功能测试（均在SYSV目录下使用）
在SYSV目录下，使用以下命令，其中使用了对应的编译选项进行编译
```
bash run_perf_test.sh
```
### 单独测试某一规模
若想单独测试某一规模，则使用命令：
```
bash run_perf_test_spec.sh type func N NRHS UPLO numThread
```
以RK函数的SINGLE类型，N=100，NRHS=50，UPLO=L，线程数48为例:(注意所有字母要大写)
```
bash run_perf_test_spec.sh SINGLE RK 100 50 L 48
```




