! Copyright (c) Huawei Technologies Co., Ltd. 2020-2020. All rights
! reserved.
! Description: Fortran API for vector math function in Kunpeng Math
! Library
! Author: xueyonghui
! Create: 2021-12-20

interface
    subroutine vsadd(len, src1, src2, dst) bind(C, name='vsadd')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src1
        real(C_FLOAT), dimension(*), intent(in) :: src2
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsadd

    subroutine vdadd(len, src1, src2, dst) bind(C, name='vdadd')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src1
        real(C_DOUBLE), dimension(*), intent(in) :: src2
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdadd

    subroutine vssub(len, src1, src2, dst) bind(C, name='vssub')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src1
        real(C_FLOAT), dimension(*), intent(in) :: src2
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vssub

    subroutine vdsub(len, src1, src2, dst) bind(C, name='vdsub')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src1
        real(C_DOUBLE), dimension(*), intent(in) :: src2
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdsub

    subroutine vsmul(len, src1, src2, dst) bind(C, name='vsmul')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src1
        real(C_FLOAT), dimension(*), intent(in) :: src2
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsmul

    subroutine vdmul(len, src1, src2, dst) bind(C, name='vdmul')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src1
        real(C_DOUBLE), dimension(*), intent(in) :: src2
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdmul

    subroutine vsdiv(len, src1, src2, dst) bind(C, name='vsdiv')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src1
        real(C_FLOAT), dimension(*), intent(in) :: src2
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsdiv

    subroutine vddiv(len, src1, src2, dst) bind(C, name='vddiv')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src1
        real(C_DOUBLE), dimension(*), intent(in) :: src2
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vddiv

    subroutine vssqr(len, src, dst) bind(C, name='vssqr')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vssqr

    subroutine vdsqr(len, src, dst) bind(C, name='vdsqr')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdsqr

    subroutine vssqrt(len, src, dst) bind(C, name='vssqrt')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vssqrt

    subroutine vdsqrt(len, src, dst) bind(C, name='vdsqrt')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdsqrt

    subroutine vsexp(len, src, dst) bind(C, name='vsexp')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsexp

    subroutine vdexp(len, src, dst) bind(C, name='vdexp')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdexp

    subroutine vsln(len, src, dst) bind(C, name='vsln')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsln

    subroutine vdln(len, src, dst) bind(C, name='vdln')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdln

    subroutine vslog10(len, src, dst) bind(C, name='vslog10')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vslog10

    subroutine vdlog10(len, src, dst) bind(C, name='vdlog10')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdlog10

    subroutine vspow(len, src1, src2, dst) bind(C, name='vspow')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src1
        real(C_FLOAT), dimension(*), intent(in) :: src2
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vspow

    subroutine vdpow(len, src1, src2, dst) bind(C, name='vdpow')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src1
        real(C_DOUBLE), dimension(*), intent(in) :: src2
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdpow

    subroutine vssin(len, src, dst) bind(C, name='vssin')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vssin

    subroutine vdsin(len, src, dst) bind(C, name='vdsin')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdsin

    subroutine vscos(len, src, dst) bind(C, name='vscos')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vscos

    subroutine vdcos(len, src, dst) bind(C, name='vdcos')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdcos

    subroutine vstan(len, src, dst) bind(C, name='vstan')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vstan

    subroutine vdtan(len, src, dst) bind(C, name='vdtan')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdtan

    subroutine vsatan(len, src, dst) bind(C, name='vsatan')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsatan

    subroutine vdatan(len, src, dst) bind(C, name='vdatan')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdatan

    subroutine vsatan2(len, srcx, srcy, dst) bind(C, name='vsatan2')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: srcx
        real(C_FLOAT), dimension(*), intent(in) :: srcy
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsatan2

    subroutine vdatan2(len, srcx, srcy, dst) bind(C, name='vdatan2')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: srcx
        real(C_DOUBLE), dimension(*), intent(in) :: srcy
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdatan2

    subroutine vssincos(len, src, sindst, cosdst) bind(C,name='vssincos')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: sindst
        real(C_FLOAT), dimension(*), intent(out) :: cosdst
    end subroutine vssincos

    subroutine vdsincos(len, src, sindst, cosdst) bind(C,name='vdsincos')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: sindst
        real(C_DOUBLE), dimension(*), intent(out) :: cosdst
    end subroutine vdsincos

    subroutine vscosh(len, src, dst) bind(C, name='vscosh')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vscosh

    subroutine vdcosh(len, src, dst) bind(C, name='vdcosh')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdcosh

    subroutine vssinh(len, src, dst) bind(C, name='vssinh')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vssinh

    subroutine vdsinh(len, src, dst) bind(C, name='vdsinh')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdsinh

    subroutine vsatanh(len, src, dst) bind(C, name='vsatanh')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsatanh

    subroutine vdatanh(len, src, dst) bind(C, name='vdatanh')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdatanh

    subroutine vsacosh(len, src, dst) bind(C, name='vsacosh')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsacosh

    subroutine vdacosh(len, src, dst) bind(C, name='vdacosh')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdacosh

    subroutine vsasinh(len, src, dst) bind(C, name='vsasinh')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vsasinh

    subroutine vdasinh(len, src, dst) bind(C, name='vdasinh')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdasinh

    subroutine vstanh(len, src, dst) bind(C, name='vstanh')
        import
        integer(C_INT), value :: len
        real(C_FLOAT), dimension(*), intent(in) :: src
        real(C_FLOAT), dimension(*), intent(out) :: dst
    end subroutine vstanh

    subroutine vdtanh(len, src, dst) bind(C, name='vdtanh')
        import
        integer(C_INT), value :: len
        real(C_DOUBLE), dimension(*), intent(in) :: src
        real(C_DOUBLE), dimension(*), intent(out) :: dst
    end subroutine vdtanh
end interface
