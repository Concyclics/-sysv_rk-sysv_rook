! Copyright (c) Huawei Technologies Co., Ltd. 2021-2021. All rights reserved.
! Description: Fortran API for Basic math function in Kunpeng Math Library
! Author: xueyonghui
! Create: 2021-05-19

interface
    real(C_FLOAT) function expf(input) bind(C, name='expf')
        import
        real(C_FLOAT), value :: input
    end function expf

    real(C_DOUBLE) function exp(input) bind(C, name='exp')
        import
        real(C_DOUBLE), value :: input
    end function exp

    real(C_FLOAT) function exp2f(input) bind(C, name='exp2f')
        import
        real(C_FLOAT), value :: input
    end function exp2f

    real(C_DOUBLE) function exp2(input) bind(C, name='exp2')
        import
        real(C_DOUBLE), value :: input
    end function exp2

    real(C_FLOAT) function logf(input) bind(C, name='logf')
        import
        real(C_FLOAT), value :: input
    end function logf

    real(C_DOUBLE) function log(input) bind(C, name='log')
        import
        real(C_DOUBLE), value :: input
    end function log

    real(C_FLOAT) function log2f(input) bind(C, name='log2f')
        import
        real(C_FLOAT), value :: input
    end function log2f

    real(C_DOUBLE) function log2(input) bind(C, name='log2')
        import
        real(C_DOUBLE), value :: input
    end function log2

    real(C_FLOAT) function log10f(input) bind(C, name='log10f')
        import
        real(C_FLOAT), value :: input
    end function log10f

    real(C_DOUBLE) function log10(input) bind(C, name='log10')
        import
        real(C_DOUBLE), value :: input
    end function log10

    real(C_FLOAT) function powf(input1, input2) bind(C, name='powf')
        import
        real(C_FLOAT), value :: input1
        real(C_FLOAT), value :: input2
    end function powf

    real(C_DOUBLE) function pow(input1, input2) bind(C, name='pow')
        import
        real(C_DOUBLE), value :: input1
        real(C_DOUBLE), value :: input2
    end function pow

    real(C_FLOAT) function sinf(input) bind(C, name='sinf')
        import
        real(C_FLOAT), value :: input
    end function sinf

    real(C_DOUBLE) function sin(input) bind(C, name='sin')
        import
        real(C_DOUBLE), value :: input
    end function sin

    real(C_FLOAT) function cosf(input) bind(C, name='cosf')
        import
        real(C_FLOAT), value :: input
    end function cosf

    real(C_DOUBLE) function cos(input) bind(C, name='cos')
        import
        real(C_DOUBLE), value :: input
    end function cos

    real(C_FLOAT) function tanf(input) bind(C, name='tanf')
        import
        real(C_FLOAT), value :: input
    end function tanf

    real(C_DOUBLE) function tan(input) bind(C, name='tan')
        import
        real(C_DOUBLE), value :: input
    end function tan

    real(C_FLOAT) function atanf(input) bind(C, name='atanf')
        import
        real(C_FLOAT), value :: input
    end function atanf

    real(C_DOUBLE) function atan(input) bind(C, name='atan')
        import
        real(C_DOUBLE), value :: input
    end function atan

    real(C_FLOAT) function atan2f(input1, input2) bind(C, name='atan2f')
        import
        real(C_FLOAT), value :: input1
        real(C_FLOAT), value :: input2
    end function atan2f

    real(C_DOUBLE) function atan2(input1, input2) bind(C, name='atan2')
        import
        real(C_DOUBLE), value :: input1
        real(C_DOUBLE), value :: input2
    end function atan2

    real(C_FLOAT) function sqrtf(input) bind(C, name='sqrtf')
        import
        real(C_FLOAT), value :: input
    end function sqrtf

    real(C_DOUBLE) function sqrt(input) bind(C, name='sqrt')
        import
        real(C_DOUBLE), value :: input
    end function sqrt

    subroutine sincosf(input, siny, cosy) bind(C, name='sincosf')
        import
        real(C_FLOAT), value :: input
        type(C_PTR), value :: siny
        type(C_PTR), value :: cosy
    end subroutine sincosf

    subroutine sincos(input, siny, cosy) bind(C, name='sincos')
        import
        real(C_DOUBLE), value :: input
        type(C_PTR), value :: siny
        type(C_PTR), value :: cosy
    end subroutine sincos
    
    real(C_FLOAT) function asinhf(input) bind(C, name='asinhf')
        import
        real(C_FLOAT), value :: input
    end function asinhf

    real(C_DOUBLE) function asinh(input) bind(C, name='asinh')
        import
        real(C_DOUBLE), value :: input
    end function asinh

    real(C_FLOAT) function acoshf(input) bind(C, name='acoshf')
        import
        real(C_FLOAT), value :: input
    end function acoshf

    real(C_DOUBLE) function acosh(input) bind(C, name='acosh')
        import
        real(C_DOUBLE), value :: input
    end function acosh

    real(C_FLOAT) function tanhf(input) bind(C, name='tanhf')
        import
        real(C_FLOAT), value :: input
    end function tanhf

    real(C_DOUBLE) function tanh(input) bind(C, name='tanh')
        import
        real(C_DOUBLE), value :: input
    end function tanh

    real(C_FLOAT) function cbrtf(input) bind(C, name='cbrtf')
        import
        real(C_FLOAT), value :: input
    end function cbrtf

    real(C_DOUBLE) function cbrt(input) bind(C, name='cbrt')
        import
        real(C_DOUBLE), value :: input
    end function cbrt

    real(C_FLOAT) function atanhf(input) bind(C, name='atanhf')
        import
        real(C_FLOAT), value :: input
    end function atanhf

    real(C_DOUBLE) function atanh(input) bind(C, name='atanh')
        import
        real(C_DOUBLE), value :: input
    end function atanh

    real(C_FLOAT) function coshf(input) bind(C, name='coshf')
        import
        real(C_FLOAT), value :: input
    end function coshf

    real(C_DOUBLE) function cosh(input) bind(C, name='cosh')
        import
        real(C_DOUBLE), value :: input
    end function cosh

    real(C_FLOAT) function sinhf(input) bind(C, name='sinhf')
        import
        real(C_FLOAT), value :: input
    end function sinhf

    real(C_DOUBLE) function sinh(input) bind(C, name='sinh')
        import
        real(C_DOUBLE), value :: input
    end function sinh

    complex(C_FLOAT_COMPLEX) function csinf(input) bind(C, name='csinf')
        import
        complex(C_FLOAT_COMPLEX), value :: input
    end function csinf

    complex(C_DOUBLE_COMPLEX) function csin(input) bind(C, name='csin')
        import
        complex(C_DOUBLE_COMPLEX), value :: input
    end function csin

    complex(C_FLOAT_COMPLEX) function ccosf(input) bind(C, name='ccosf')
        import
        complex(C_FLOAT_COMPLEX), value :: input
    end function ccosf

    complex(C_DOUBLE_COMPLEX) function ccos(input) bind(C, name='ccos')
        import
        complex(C_DOUBLE_COMPLEX), value :: input
    end function ccos

    complex(C_FLOAT_COMPLEX) function ctanf(input) bind(C, name='ctanf')
        import
        complex(C_FLOAT_COMPLEX), value :: input
    end function ctanf

    complex(C_DOUBLE_COMPLEX) function ctan(input) bind(C, name='ctan')
        import
        complex(C_DOUBLE_COMPLEX), value :: input
    end function ctan

    complex(C_FLOAT_COMPLEX) function catanf(input) bind(C, name='catanf')
        import
        complex(C_FLOAT_COMPLEX), value :: input
    end function catanf

    complex(C_DOUBLE_COMPLEX) function catan(input) bind(C, name='catan')
        import
        complex(C_DOUBLE_COMPLEX), value :: input
    end function catan

end interface