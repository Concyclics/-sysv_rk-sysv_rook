! Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights reserved.
! Description: Fortran API for Kunpeng Math Library FFT
! Author: xueyonghui
! Create: 2021-01-08

integer(C_INT), parameter :: KML_FFT_FORWARD = -1
integer(C_INT), parameter :: KML_FFT_BACKWARD = +1
integer(C_INT), parameter :: KML_FFT_MEASURE = 0
integer(C_INT), parameter :: KML_FFT_PATIENT = 32
integer(C_INT), parameter :: KML_FFT_ESTIMATE = 64
integer(C_INT), parameter :: KML_FFT_R2HC = 0
integer(C_INT), parameter :: KML_FFT_HC2R = 1
integer(C_INT), parameter :: KML_FFT_DHT = 2
integer(C_INT), parameter :: KML_FFT_REDFT00 = 3
integer(C_INT), parameter :: KML_FFT_REDFT01 = 4
integer(C_INT), parameter :: KML_FFT_REDFT10 = 5
integer(C_INT), parameter :: KML_FFT_REDFT11 = 6
integer(C_INT), parameter :: KML_FFT_RODFT00 = 7
integer(C_INT), parameter :: KML_FFT_RODFT01 = 8
integer(C_INT), parameter :: KML_FFT_RODFT10 = 9
integer(C_INT), parameter :: KML_FFT_RODFT11 = 10

type, bind(C) :: kml_fft_iodim
    integer(C_INT) n
    integer(C_INT) is
    integer(C_INT) os
end type kml_fft_iodim

type, bind(C) :: kml_fft_iodim64
    integer(C_INTPTR_T) n
    integer(C_INTPTR_T) is
    integer(C_INTPTR_T) os
end type kml_fft_iodim64

type, bind(C) :: kml_fft_complex
    real(C_DOUBLE) r
    real(C_DOUBLE) i
end type kml_fft_complex


interface
    type(C_PTR) function kml_fft_plan_dft_1d(n, in, out, sign, flags) bind(C, name='kml_fft_plan_dft_1d')
        import
        integer(C_INT), value :: n
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_1d

    type(C_PTR) function kml_fft_plan_dft_2d(n0, n1, in, out, sign, flags) bind(C, name='kml_fft_plan_dft_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_2d

    type(C_PTR) function kml_fft_plan_dft_3d(n0, n1, n2, in, out, sign, flags) bind(C, name='kml_fft_plan_dft_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_3d

    type(C_PTR) function kml_fft_plan_dft(rank, n, in, out, sign, flags) bind(C, name='kml_fft_plan_dft')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft

    type(C_PTR) function kml_fft_plan_many_dft(&
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags) &
            bind(C, name='kml_fft_plan_many_dft')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        type(kml_fft_complex), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist 
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_many_dft

    type(C_PTR) function kml_fft_plan_dft_r2c_1d(n, in, out, flags) bind(C, name='kml_fft_plan_dft_r2c_1d')
        import
        integer(C_INT), value :: n
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_r2c_1d

    type(C_PTR) function kml_fft_plan_dft_r2c_2d(n0, n1, in, out, flags) bind(C, name='kml_fft_plan_dft_r2c_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_r2c_2d

    type(C_PTR) function kml_fft_plan_dft_r2c_3d(n0, n1, n2, in, out, flags) bind(C, name='kml_fft_plan_dft_r2c_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_r2c_3d

    type(C_PTR) function kml_fft_plan_dft_r2c(rank, n, in, out, flags) bind(C, name='kml_fft_plan_dft_r2c')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_r2c

    type(C_PTR) function kml_fft_plan_many_dft_r2c( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) &
            bind(C, name='kml_fft_plan_many_dft_r2c')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        real(C_DOUBLE), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist 
        integer(C_INT), value :: flags
    end function kml_fft_plan_many_dft_r2c

    type(C_PTR) function kml_fft_plan_dft_c2r_1d(n, in, out, flags) bind(C, name='kml_fft_plan_dft_c2r_1d')
        import
        integer(C_INT), value :: n
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_c2r_1d

    type(C_PTR) function kml_fft_plan_dft_c2r_2d(n0, n1, in, out, flags) bind(C, name='kml_fft_plan_dft_c2r_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_c2r_2d

    type(C_PTR) function kml_fft_plan_dft_c2r_3d(n0, n1, n2, in, out, flags) bind(C, name='kml_fft_plan_dft_c2r_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_c2r_3d

    type(C_PTR) function kml_fft_plan_dft_c2r(rank, n, in, out, flags) bind(C, name='kml_fft_plan_dft_c2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_dft_c2r

    type(C_PTR) function kml_fft_plan_many_dft_c2r( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) &
            bind(C, name='kml_fft_plan_many_dft_c2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        type(kml_fft_complex), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist 
        integer(C_INT), value :: flags
    end function kml_fft_plan_many_dft_c2r

    type(C_PTR) function kml_fft_plan_r2r_1d(n, in, out, kind, flags) bind(C, name='kml_fft_plan_r2r_1d')
        import
        integer(C_INT), value :: n
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: kind
        integer(C_INT), value :: flags
    end function kml_fft_plan_r2r_1d

    type(C_PTR) function kml_fft_plan_r2r_2d(n0, n1, in, out, kind0, kind1, flags) &
            bind(C, name='kml_fft_plan_r2r_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: kind0
        integer(C_INT), value :: kind1
        integer(C_INT), value :: flags
    end function kml_fft_plan_r2r_2d

    type(C_PTR) function kml_fft_plan_r2r_3d(n0, n1, n2, in, out, kind0, kind1, kind2, flags) & 
            bind(C, name='kml_fft_plan_r2r_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: kind0
        integer(C_INT), value :: kind1
        integer(C_INT), value :: kind2
        integer(C_INT), value :: flags
    end function kml_fft_plan_r2r_3d

    type(C_PTR) function kml_fft_plan_r2r(rank, n, in, out, kind, flags) bind(C, name='kml_fft_plan_r2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: kind
        integer(C_INT), value :: flags
    end function kml_fft_plan_r2r

    type(C_PTR) function kml_fft_plan_many_r2r( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, kind, flags) &
            bind(C, name='kml_fft_plan_many_r2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        real(C_DOUBLE), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist
        integer(C_INT), dimension(*), intent(in) :: kind 
        integer(C_INT), value :: flags
    end function kml_fft_plan_many_r2r

    type(C_PTR) function kml_fft_plan_guru_dft(rank, dims, howmany_rank, howmany_dims, in, out, sign, flags) &
            bind(C, name='kml_fft_plan_guru_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_dft

    type(C_PTR) function kml_fft_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, ri, ii, ro, io, flags) &
            bind(C, name='kml_fft_plan_guru_split_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: ri
        real(C_DOUBLE), dimension(*), intent(out) :: ii
        real(C_DOUBLE), dimension(*), intent(out) :: ro
        real(C_DOUBLE), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_split_dft

    type(C_PTR) function kml_fft_plan_guru64_dft(rank, dims, howmany_rank, howmany_dims, in, out, sign, flags) &
            bind(C, name='kml_fft_plan_guru64_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_dft

    type(C_PTR) function kml_fft_plan_guru64_split_dft(rank, dims, howmany_rank, howmany_dims, ri, ii, ro, io, flags) &
            bind(C, name='kml_fft_plan_guru64_split_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: ri
        real(C_DOUBLE), dimension(*), intent(out) :: ii
        real(C_DOUBLE), dimension(*), intent(out) :: ro
        real(C_DOUBLE), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_split_dft

    type(C_PTR) function kml_fft_plan_guru_dft_r2c(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fft_plan_guru_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_dft_r2c

    type(C_PTR) function kml_fft_plan_guru_split_dft_r2c( &
        rank, dims, howmany_rank, howmany_dims, in, ro, io, flags) &
            bind(C, name='kml_fft_plan_guru_split_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: ro
        real(C_DOUBLE), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_split_dft_r2c

    type(C_PTR) function kml_fft_plan_guru64_dft_r2c(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fft_plan_guru64_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_dft_r2c

    type(C_PTR) function kml_fft_plan_guru64_split_dft_r2c( &
        rank, dims, howmany_rank, howmany_dims, in, ro, io, flags) bind(C, name='kml_fft_plan_guru64_split_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: ro
        real(C_DOUBLE), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_split_dft_r2c

    type(C_PTR) function kml_fft_plan_guru_dft_c2r(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fft_plan_guru_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_dft_c2r

    type(C_PTR) function kml_fft_plan_guru_split_dft_c2r(rank, dims, howmany_rank, howmany_dims, ri, ii, out, flags) &
            bind(C, name='kml_fft_plan_guru_split_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: ri
        real(C_DOUBLE), dimension(*), intent(out) :: ii
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_split_dft_c2r

    type(C_PTR) function kml_fft_plan_guru64_dft_c2r(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fft_plan_guru64_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_dft_c2r

    type(C_PTR) function kml_fft_plan_guru64_split_dft_c2r( &
        rank, dims, howmany_rank, howmany_dims, ri, ii, out, flags) &
            bind(C, name='kml_fft_plan_guru64_split_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: ri
        real(C_DOUBLE), dimension(*), intent(out) :: ii
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_split_dft_c2r

    type(C_PTR) function kml_fft_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, in, out, kind, flags) &
            bind(C, name='kml_fft_plan_guru_r2r')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: kind
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru_r2r

    type(C_PTR) function kml_fft_plan_guru64_r2r(rank, dims, howmany_rank, howmany_dims, in, out, kind, flags) &
            bind(C, name='kml_fft_plan_guru64_r2r')
        import
        integer(C_INT), value :: rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fft_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: kind
        integer(C_INT), value :: flags
    end function kml_fft_plan_guru64_r2r

    subroutine kml_fft_execute(p) bind(C, name='kml_fft_execute')
        import
        type(C_PTR), value :: p
    end subroutine kml_fft_execute

    subroutine kml_fft_execute_dft(p, in, out) bind(C, name='kml_fft_execute_dft')
        import
        type(C_PTR), value :: p
        type(kml_fft_complex), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
    end subroutine kml_fft_execute_dft

    subroutine kml_fft_execute_dft_r2c(p, in, out) bind(C, name='kml_fft_execute_dft_r2c')
        import
        type(C_PTR), value :: p
        real(C_DOUBLE), dimension(*), intent(out) :: in
        type(kml_fft_complex), dimension(*), intent(out) :: out
    end subroutine kml_fft_execute_dft_r2c

    subroutine kml_fft_execute_dft_c2r(p, in, out) bind(C, name='kml_fft_execute_dft_c2r')
        import
        type(C_PTR), value :: p
        type(kml_fft_complex), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine kml_fft_execute_dft_c2r

    subroutine kml_fft_execute_split_dft(p, ri, ii, ro, io) bind(C, name='kml_fft_execute_split_dft')
        import
        type(C_PTR), value :: p
        real(C_DOUBLE), dimension(*), intent(out) :: ri
        real(C_DOUBLE), dimension(*), intent(out) :: ii
        real(C_DOUBLE), dimension(*), intent(out) :: ro
        real(C_DOUBLE), dimension(*), intent(out) :: io
    end subroutine kml_fft_execute_split_dft

    subroutine kml_fft_execute_split_dft_r2c(p, in, ro, io) bind(C, name='kml_fft_execute_split_dft_r2c')
        import
        type(C_PTR), value :: p
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: ro
        real(C_DOUBLE), dimension(*), intent(out) :: io
    end subroutine kml_fft_execute_split_dft_r2c

    subroutine kml_fft_execute_split_dft_c2r(p, ri, ii, out) bind(C, name='kml_fft_execute_split_dft_c2r')
        import
        type(C_PTR), value :: p
        real(C_DOUBLE), dimension(*), intent(out) :: ri
        real(C_DOUBLE), dimension(*), intent(out) :: ii
        real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine kml_fft_execute_split_dft_c2r

    subroutine kml_fft_execute_r2r(p, in, out) bind(C, name='kml_fft_execute_r2r')
        import
        type(C_PTR), value :: p
        real(C_DOUBLE), dimension(*), intent(out) :: in
        real(C_DOUBLE), dimension(*), intent(out) :: out
    end subroutine kml_fft_execute_r2r

    subroutine kml_fft_destroy_plan(p) bind(C, name='kml_fft_destroy_plan')
        import
        type(C_PTR), value :: p
    end subroutine kml_fft_destroy_plan

    subroutine kml_fft_plan_with_nthreads(nthreads) bind(C, name='kml_fft_plan_with_nthreads')
        import
        integer(C_INT), value :: nthreads
    end subroutine kml_fft_plan_with_nthreads

    integer(C_INT) function kml_fft_init_threads() bind(C, name='kml_fft_init_threads')
        import
    end function kml_fft_init_threads

    subroutine kml_fft_cleanup_threads() bind(C, name='kml_fft_cleanup_threads')
        import
    end subroutine kml_fft_cleanup_threads

    type(C_PTR) function kml_fft_malloc(n) bind(C, name='kml_fft_malloc')
        import
        integer(C_SIZE_T), value :: n
    end function kml_fft_malloc

    subroutine kml_fft_free(p) bind(C, name='kml_fft_free')
        import
        type(C_PTR), value :: p
    end subroutine kml_fft_free
end interface

type, bind(C) :: kml_fftf_iodim
    integer(C_INT) n
    integer(C_INT) is
    integer(C_INT) os
end type kml_fftf_iodim

type, bind(C) :: kml_fftf_iodim64
    integer(C_INTPTR_T) n
    integer(C_INTPTR_T) is
    integer(C_INTPTR_T) os
end type kml_fftf_iodim64

type, bind(C) :: kml_fftf_complex
    real(C_FLOAT) r
    real(C_FLOAT) i
end type kml_fftf_complex

interface
    type(C_PTR) function kml_fftf_plan_dft_1d(n, in, out, sign, flags) bind(C, name='kml_fftf_plan_dft_1d')
        import
        integer(C_INT), value :: n
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_1d

    type(C_PTR) function kml_fftf_plan_dft_2d(n0, n1, in, out, sign, flags) bind(C, name='kml_fftf_plan_dft_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_2d

    type(C_PTR) function kml_fftf_plan_dft_3d(n0, n1, n2, in, out, sign, flags) bind(C, name='kml_fftf_plan_dft_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_3d

    type(C_PTR) function kml_fftf_plan_dft(rank, n, in, out, sign, flags) bind(C, name='kml_fftf_plan_dft')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft

    type(C_PTR) function kml_fftf_plan_many_dft( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags) &
            bind(C, name='kml_fftf_plan_many_dft')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist 
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_many_dft

    type(C_PTR) function kml_fftf_plan_dft_r2c_1d(n, in, out, flags) bind(C, name='kml_fftf_plan_dft_r2c_1d')
        import
        integer(C_INT), value :: n
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_r2c_1d

    type(C_PTR) function kml_fftf_plan_dft_r2c_2d(n0, n1, in, out, flags) bind(C, name='kml_fftf_plan_dft_r2c_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_r2c_2d

    type(C_PTR) function kml_fftf_plan_dft_r2c_3d(n0, n1, n2, in, out, flags) bind(C, name='kml_fftf_plan_dft_r2c_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_r2c_3d

    type(C_PTR) function kml_fftf_plan_dft_r2c(rank, n, in, out, flags) bind(C, name='kml_fftf_plan_dft_r2c')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_r2c

    type(C_PTR) function kml_fftf_plan_many_dft_r2c( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) &
            bind(C, name='kml_fftf_plan_many_dft_r2c')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        real(C_FLOAT), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist 
        integer(C_INT), value :: flags
    end function kml_fftf_plan_many_dft_r2c

    type(C_PTR) function kml_fftf_plan_dft_c2r_1d(n, in, out, flags) bind(C, name='kml_fftf_plan_dft_c2r_1d')
        import
        integer(C_INT), value :: n
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_c2r_1d

    type(C_PTR) function kml_fftf_plan_dft_c2r_2d(n0, n1, in, out, flags) bind(C, name='kml_fftf_plan_dft_c2r_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_c2r_2d

    type(C_PTR) function kml_fftf_plan_dft_c2r_3d(n0, n1, n2, in, out, flags) bind(C, name='kml_fftf_plan_dft_c2r_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_c2r_3d

    type(C_PTR) function kml_fftf_plan_dft_c2r(rank, n, in, out, flags) bind(C, name='kml_fftf_plan_dft_c2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_dft_c2r

    type(C_PTR) function kml_fftf_plan_many_dft_c2r( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, flags) &
            bind(C, name='kml_fftf_plan_many_dft_c2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist 
        integer(C_INT), value :: flags
    end function kml_fftf_plan_many_dft_c2r

    type(C_PTR) function kml_fftf_plan_r2r_1d(n, in, out, kind, flags) bind(C, name='kml_fftf_plan_r2r_1d')
        import
        integer(C_INT), value :: n
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: kind
        integer(C_INT), value :: flags
    end function kml_fftf_plan_r2r_1d

    type(C_PTR) function kml_fftf_plan_r2r_2d(n0, n1, in, out, kind0, kind1, flags) &
            bind(C, name='kml_fftf_plan_r2r_2d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: kind0
        integer(C_INT), value :: kind1
        integer(C_INT), value :: flags
    end function kml_fftf_plan_r2r_2d

    type(C_PTR) function kml_fftf_plan_r2r_3d(n0, n1, n2, in, out, kind0, kind1, kind2, flags) & 
            bind(C, name='kml_fftf_plan_r2r_3d')
        import
        integer(C_INT), value :: n0
        integer(C_INT), value :: n1
        integer(C_INT), value :: n2
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: kind0
        integer(C_INT), value :: kind1
        integer(C_INT), value :: kind2
        integer(C_INT), value :: flags
    end function kml_fftf_plan_r2r_3d

    type(C_PTR) function kml_fftf_plan_r2r(rank, n, in, out, kind, flags) bind(C, name='kml_fftf_plan_r2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: kind
        integer(C_INT), value :: flags
    end function kml_fftf_plan_r2r

    type(C_PTR) function kml_fftf_plan_many_r2r( &
        rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, kind, flags) &
            bind(C, name='kml_fftf_plan_many_r2r')
        import
        integer(C_INT), value :: rank
        integer(C_INT), dimension(*), intent(in) :: n
        integer(C_INT), value :: howmany
        real(C_FLOAT), dimension(*), intent(out) :: in
        integer(C_INT), dimension(*), intent(in) :: inembed
        integer(C_INT), value :: istride
        integer(C_INT), value :: idist 
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: onembed
        integer(C_INT), value :: ostride
        integer(C_INT), value :: odist
        integer(C_INT), dimension(*), intent(in) :: kind 
        integer(C_INT), value :: flags
    end function kml_fftf_plan_many_r2r

    type(C_PTR) function kml_fftf_plan_guru_dft(rank, dims, howmany_rank, howmany_dims, in, out, sign, flags) &
            bind(C, name='kml_fftf_plan_guru_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_dft

    type(C_PTR) function kml_fftf_plan_guru_split_dft(rank, dims, howmany_rank, howmany_dims, ri, ii, ro, io, flags) &
            bind(C, name='kml_fftf_plan_guru_split_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: ri
        real(C_FLOAT), dimension(*), intent(out) :: ii
        real(C_FLOAT), dimension(*), intent(out) :: ro
        real(C_FLOAT), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_split_dft

    type(C_PTR) function kml_fftf_plan_guru64_dft(rank, dims, howmany_rank, howmany_dims, in, out, sign, flags) &
            bind(C, name='kml_fftf_plan_guru64_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: sign
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_dft

    type(C_PTR) function kml_fftf_plan_guru64_split_dft(rank, dims, howmany_rank, howmany_dims, ri, ii, ro, io, flags) &
            bind(C, name='kml_fftf_plan_guru64_split_dft')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: ri
        real(C_FLOAT), dimension(*), intent(out) :: ii
        real(C_FLOAT), dimension(*), intent(out) :: ro
        real(C_FLOAT), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_split_dft

    type(C_PTR) function kml_fftf_plan_guru_dft_r2c(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fftf_plan_guru_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_dft_r2c

    type(C_PTR) function kml_fftf_plan_guru_split_dft_r2c( &
        rank, dims, howmany_rank, howmany_dims, in, ro, io, flags) &
            bind(C, name='kml_fftf_plan_guru_split_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: ro
        real(C_FLOAT), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_split_dft_r2c

    type(C_PTR) function kml_fftf_plan_guru64_dft_r2c(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fftf_plan_guru64_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_dft_r2c

    type(C_PTR) function kml_fftf_plan_guru64_split_dft_r2c( &
        rank, dims, howmany_rank, howmany_dims, in, ro, io, flags) bind(C, name='kml_fftf_plan_guru64_split_dft_r2c')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: ro
        real(C_FLOAT), dimension(*), intent(out) :: io
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_split_dft_r2c

    type(C_PTR) function kml_fftf_plan_guru_dft_c2r(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fftf_plan_guru_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_dft_c2r

    type(C_PTR) function kml_fftf_plan_guru_split_dft_c2r(rank, dims, howmany_rank, howmany_dims, ri, ii, out, flags) &
            bind(C, name='kml_fftf_plan_guru_split_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: ri
        real(C_FLOAT), dimension(*), intent(out) :: ii
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_split_dft_c2r

    type(C_PTR) function kml_fftf_plan_guru64_dft_c2r(rank, dims, howmany_rank, howmany_dims, in, out, flags) &
            bind(C, name='kml_fftf_plan_guru64_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_dft_c2r

    type(C_PTR) function kml_fftf_plan_guru64_split_dft_c2r( &
        rank, dims, howmany_rank, howmany_dims, ri, ii, out, flags) &
            bind(C, name='kml_fftf_plan_guru64_split_dft_c2r')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: ri
        real(C_FLOAT), dimension(*), intent(out) :: ii
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_split_dft_c2r

    type(C_PTR) function kml_fftf_plan_guru_r2r(rank, dims, howmany_rank, howmany_dims, in, out, kind, flags) &
            bind(C, name='kml_fftf_plan_guru_r2r')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: kind
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru_r2r

    type(C_PTR) function kml_fftf_plan_guru64_r2r(rank, dims, howmany_rank, howmany_dims, in, out, kind, flags) &
            bind(C, name='kml_fftf_plan_guru64_r2r')
        import
        integer(C_INT), value :: rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: dims
        integer(C_INT), value :: howmany_rank
        type(kml_fftf_iodim64), dimension(*), intent(in) :: howmany_dims
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
        integer(C_INT), dimension(*), intent(in) :: kind
        integer(C_INT), value :: flags
    end function kml_fftf_plan_guru64_r2r

    subroutine kml_fftf_execute(p) bind(C, name='kml_fftf_execute')
        import
        type(C_PTR), value :: p
    end subroutine kml_fftf_execute

    subroutine kml_fftf_execute_dft(p, in, out) bind(C, name='kml_fftf_execute_dft')
        import
        type(C_PTR), value :: p
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
    end subroutine kml_fftf_execute_dft

    subroutine kml_fftf_execute_dft_r2c(p, in, out) bind(C, name='kml_fftf_execute_dft_r2c')
        import
        type(C_PTR), value :: p
        real(C_FLOAT), dimension(*), intent(out) :: in
        type(kml_fftf_complex), dimension(*), intent(out) :: out
    end subroutine kml_fftf_execute_dft_r2c

    subroutine kml_fftf_execute_dft_c2r(p, in, out) bind(C, name='kml_fftf_execute_dft_c2r')
        import
        type(C_PTR), value :: p
        type(kml_fftf_complex), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine kml_fftf_execute_dft_c2r

    subroutine kml_fftf_execute_split_dft(p, ri, ii, ro, io) bind(C, name='kml_fftf_execute_split_dft')
        import
        type(C_PTR), value :: p
        real(C_FLOAT), dimension(*), intent(out) :: ri
        real(C_FLOAT), dimension(*), intent(out) :: ii
        real(C_FLOAT), dimension(*), intent(out) :: ro
        real(C_FLOAT), dimension(*), intent(out) :: io
    end subroutine kml_fftf_execute_split_dft

    subroutine kml_fftf_execute_split_dft_r2c(p, in, ro, io) bind(C, name='kml_fftf_execute_split_dft_r2c')
        import
        type(C_PTR), value :: p
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: ro
        real(C_FLOAT), dimension(*), intent(out) :: io
    end subroutine kml_fftf_execute_split_dft_r2c

    subroutine kml_fftf_execute_split_dft_c2r(p, ri, ii, out) bind(C, name='kml_fftf_execute_split_dft_c2r')
        import
        type(C_PTR), value :: p
        real(C_FLOAT), dimension(*), intent(out) :: ri
        real(C_FLOAT), dimension(*), intent(out) :: ii
        real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine kml_fftf_execute_split_dft_c2r

    subroutine kml_fftf_execute_r2r(p, in, out) bind(C, name='kml_fftf_execute_r2r')
        import
        type(C_PTR), value :: p
        real(C_FLOAT), dimension(*), intent(out) :: in
        real(C_FLOAT), dimension(*), intent(out) :: out
    end subroutine kml_fftf_execute_r2r

    subroutine kml_fftf_destroy_plan(p) bind(C, name='kml_fftf_destroy_plan')
        import
        type(C_PTR), value :: p
    end subroutine kml_fftf_destroy_plan

    subroutine kml_fftf_plan_with_nthreads(nthreads) bind(C, name='kml_fftf_plan_with_nthreads')
        import
        integer(C_INT), value :: nthreads
    end subroutine kml_fftf_plan_with_nthreads

    integer(C_INT) function kml_fftf_init_threads() bind(C, name='kml_fftf_init_threads')
        import
    end function kml_fftf_init_threads

    subroutine kml_fftf_cleanup_threads() bind(C, name='kml_fftf_cleanup_threads')
        import
    end subroutine kml_fftf_cleanup_threads

    type(C_PTR) function kml_fftf_malloc(n) bind(C, name='kml_fftf_malloc')
        import
        integer(C_SIZE_T), value :: n
    end function kml_fftf_malloc

    subroutine kml_fftf_free(p) bind(C, name='kml_fftf_free')
        import
        type(C_PTR), value :: p
    end subroutine kml_fftf_free
end interface

