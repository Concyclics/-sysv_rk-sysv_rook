! Copyright (c) Huawei Technologies Co., Ltd. 2020-2021. All rights reserved.
! Description: Fortran API for Kunpeng Math Library Sparse BLAS
! Author: xueyonghui
! Create: 2021-01-09

integer(C_INT), parameter :: KML_SPARSE_STATUS_SUCCESS = 0
integer(C_INT), parameter :: KML_SPARSE_STATUS_NOT_INITIALIZED = 1
integer(C_INT), parameter :: KML_SPARSE_STATUS_ALLOC_FAILED = 2
integer(C_INT), parameter :: KML_SPARSE_STATUS_INVALID_VALUE = 3
integer(C_INT), parameter :: KML_SPARSE_STATUS_EXECUTION_FAILED = 4
integer(C_INT), parameter :: KML_SPARSE_STATUS_INTERNAL_ERROR = 5
integer(C_INT), parameter :: KML_SPARSE_STATUS_NOT_SUPPORTED = 6

integer(C_INT), parameter :: KML_SPARSE_OPERATION_NON_TRANSPOSE = 0
integer(C_INT), parameter :: KML_SPARSE_OPERATION_TRANSPOSE = 1
integer(C_INT), parameter :: KML_SPARSE_OPERATION_CONJUGATE_TRANSPOSE = 2
integer(C_INT), parameter :: KML_SPARSE_OPERATION_NUM = 3

integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_GENERAL = 0
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_SYMMETRIC = 1
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_HERMITIAN = 2
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_TRIANGULAR = 3
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_DIAGONAL = 4
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR = 5
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_BLOCK_DIAGONAL = 6
integer(C_INT), parameter :: KML_SPARSE_MATRIX_TYPE_NUM = 7

integer(C_INT), parameter :: KML_SPARSE_INDEX_BASE_ZERO = 0
integer(C_INT), parameter :: KML_SPARSE_INDEX_BASE_ONE = 1

integer(C_INT), parameter :: KML_SPARSE_FILL_MODE_LOWER = 0
integer(C_INT), parameter :: KML_SPARSE_FILL_MODE_UPPER = 1

integer(C_INT), parameter :: KML_SPARSE_DIAG_NON_UNIT = 0
integer(C_INT), parameter :: KML_SPARSE_DIAG_UNIT = 1

integer(C_INT), parameter :: KML_SPARSE_LAYOUT_ROW_MAJOR = 0
integer(C_INT), parameter :: KML_SPARSE_LAYOUT_COLUMN_MAJOR = 1

integer(C_INT), parameter :: KML_SPARSE_FORMAT_COO = 0
integer(C_INT), parameter :: KML_SPARSE_FORMAT_CSR = 1
integer(C_INT), parameter :: KML_SPARSE_FORMAT_CSC = 2
integer(C_INT), parameter :: KML_SPARSE_FORMAT_BSR = 3
integer(C_INT), parameter :: KML_SPARSE_FORMAT_SKY = 4
integer(C_INT), parameter :: KML_SPARSE_FORMAT_DIA = 5

integer(C_INT), parameter :: KML_SPARSE_DATATYPE_FLOAT = 0
integer(C_INT), parameter :: KML_SPARSE_DATATYPE_DOUBLE = 1
integer(C_INT), parameter :: KML_SPARSE_DATATYPE_FLOAT_COMPLEX = 2
integer(C_INT), parameter :: KML_SPARSE_DATATYPE_DOUBLE_COMPLEX = 3

type, bind(C) :: kml_sparse_matrix_t
    type(C_PTR) mat
    integer(C_INT) format
    integer(C_INT) datatype
end type kml_sparse_matrix_t

type, bind(C) :: kml_matrix_descr
    integer(C_INT) type
    integer(C_INT) mode
    integer(C_INT) diag
end type kml_matrix_descr

type, bind(C) :: KML_Complex8
    real(C_FLOAT) real
    real(C_FLOAT) imag
end type KML_Complex8

type, bind(C) :: KML_Complex16
    real(C_DOUBLE) real
    real(C_DOUBLE) imag
end type KML_Complex16

interface
    integer(C_INT) function kml_sparse_scsrmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_scsrmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), value :: beta
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scsrmv

    integer(C_INT) function kml_sparse_dcsrmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_dcsrmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), value :: beta
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcsrmv

    integer(C_INT) function kml_sparse_ccsrmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_ccsrmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), value :: beta
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccsrmv

    integer(C_INT) function kml_sparse_zcsrmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_zcsrmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), value :: beta
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcsrmv

    integer(C_INT) function kml_sparse_scsrsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_scsrsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scsrsv

    integer(C_INT) function kml_sparse_dcsrsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_dcsrsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcsrsv

    integer(C_INT) function kml_sparse_ccsrsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_ccsrsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccsrsv

    integer(C_INT) function kml_sparse_zcsrsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_zcsrsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcsrsv

    integer(C_INT) function kml_sparse_scsrmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) & 
                bind(C, name='kml_sparse_scsrmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_FLOAT), value :: beta
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_scsrmm

    integer(C_INT) function kml_sparse_dcsrmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) &
                bind(C, name='kml_sparse_dcsrmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_DOUBLE), value :: beta
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_dcsrmm

    integer(C_INT) function kml_sparse_ccsrmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc)  &
                bind(C, name='kml_sparse_ccsrmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex8), value :: beta
        type(KML_Complex8), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_ccsrmm

    integer(C_INT) function kml_sparse_zcsrmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) & 
                bind(C, name='kml_sparse_zcsrmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex16), value :: beta
        type(KML_Complex16), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_zcsrmm 
    
    integer(C_INT) function kml_sparse_scsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_scsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scsrgemv

    integer(C_INT) function kml_sparse_dcsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_dcsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcsrgemv

    integer(C_INT) function kml_sparse_ccsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_ccsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccsrgemv

    integer(C_INT) function kml_sparse_zcsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_zcsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcsrgemv    

    integer(C_INT) function kml_sparse_scsrsymv(uplo, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_scsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scsrsymv

    integer(C_INT) function kml_sparse_dcsrsymv(uplo, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_dcsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcsrsymv

    integer(C_INT) function kml_sparse_ccsrsymv(uplo, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_ccsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccsrsymv

    integer(C_INT) function kml_sparse_zcsrsymv(uplo, m, a, ia, ja, x, y) &
            bind(C, name='kml_sparse_zcsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcsrsymv

    integer(C_INT) function kml_csparse_scsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_csparse_scsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_csparse_scsrgemv

    integer(C_INT) function kml_csparse_dcsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_csparse_dcsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_csparse_dcsrgemv

    integer(C_INT) function kml_csparse_ccsrgemv(opt, m, a, ia, ja, x, y) &
            bind(C, name='kml_csparse_ccsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_csparse_ccsrgemv

    integer(C_INT) function kml_csparse_zcsrgemv(opt, m, a, ia, ja, x, y) & 
            bind(C, name='kml_csparse_zcsrgemv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_csparse_zcsrgemv    

    integer(C_INT) function kml_csparse_scsrsymv(uplo, m, a, ia, ja, x, y) & 
            bind(C, name='kml_csparse_scsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_csparse_scsrsymv

    integer(C_INT) function kml_csparse_dcsrsymv(uplo, m, a, ia, ja, x, y) &
            bind(C, name='kml_csparse_dcsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_csparse_dcsrsymv

    integer(C_INT) function kml_csparse_ccsrsymv(uplo, m, a, ia, ja, x, y)  &
            bind(C, name='kml_csparse_ccsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_csparse_ccsrsymv

    integer(C_INT) function kml_csparse_zcsrsymv(uplo, m, a, ia, ja, x, y) & 
            bind(C, name='kml_csparse_zcsrsymv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: m
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_csparse_zcsrsymv 
    
    integer(C_INT) function kml_sparse_scsrmultd(opt, m, n, k, a, ja, ia, b, jb, ib, c, ldc) &
            bind(C, name='kml_sparse_scsrmultd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_scsrmultd

    integer(C_INT) function kml_sparse_dcsrmultd(opt, m, n, k, a, ja, ia, b, jb, ib, c, ldc) &
            bind(C, name='kml_sparse_dcsrmultd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_dcsrmultd

    integer(C_INT) function kml_sparse_ccsrmultd(opt, m, n, k, a, ja, ia, b, jb, ib, c, ldc) &
            bind(C, name='kml_sparse_ccsrmultd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        type(KML_Complex8) , dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_ccsrmultd

    integer(C_INT) function kml_sparse_zcsrmultd(opt, m, n, k, a, ja, ia, b, jb, ib, c, ldc) &
            bind(C, name='kml_sparse_zcsrmultd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        type(KML_Complex16) , dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_zcsrmultd
    
    integer(C_INT) function kml_sparse_scsrmultcsr( &
            opt, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_scsrmultcsr')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_scsrmultcsr

    integer(C_INT) function kml_sparse_dcsrmultcsr( &
            opt, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_dcsrmultcsr')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_dcsrmultcsr

    integer(C_INT) function kml_sparse_ccsrmultcsr( &
            opt, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_ccsrmultcsr')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        type(KML_Complex8), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_ccsrmultcsr

    integer(C_INT) function kml_sparse_zcsrmultcsr( &
            opt, request, sort, m, n, k, a, ja, ia, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_zcsrmultcsr')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        type(KML_Complex16), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_zcsrmultcsr
    
    integer(C_INT) function kml_sparse_scsradd( &
            opt, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_scsradd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        real(C_FLOAT), value :: beta
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_scsradd

    integer(C_INT) function kml_sparse_dcsradd( &
            opt, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_dcsradd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        real(C_DOUBLE), value :: beta
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_dcsradd

    integer(C_INT) function kml_sparse_ccsradd( &
            opt, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_ccsradd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        type(KML_Complex8), value :: beta
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        type(KML_Complex8), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_ccsradd

    integer(C_INT) function kml_sparse_zcsradd( &
            opt, request, sort, m, n, a, ja, ia, beta, b, jb, ib, c, jc, ic, nzmax) bind(C, name='kml_sparse_zcsradd')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: request
        integer(C_INT), value :: sort
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ja
        integer(C_INT), dimension(*), intent(in) :: ia
        type(KML_Complex16), value :: beta
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), dimension(*), intent(in) :: jb
        integer(C_INT), dimension(*), intent(in) :: ib
        type(KML_Complex16), dimension(*), intent(out) :: c
        integer(C_INT), dimension(*), intent(out) :: jc
        integer(C_INT), dimension(*), intent(out) :: ic
        integer(C_INT), value :: nzmax
    end function kml_sparse_zcsradd 

    integer(C_INT) function kml_sparse_saxpyi(nz, a, x, indx, y) bind(C, name='kml_sparse_saxpyi')
        import
        integer(C_INT), value :: nz
        real(C_FLOAT), value :: a
        real(C_FLOAT), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_saxpyi    

    integer(C_INT) function kml_sparse_daxpyi(nz, a, x, indx, y) bind(C, name='kml_sparse_daxpyi')
        import
        integer(C_INT), value :: nz
        real(C_DOUBLE), value :: a
        real(C_DOUBLE), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_daxpyi    

    integer(C_INT) function kml_sparse_caxpyi(nz, a, x, indx, y) bind(C, name='kml_sparse_caxpyi')
        import
        integer(C_INT), value :: nz
        type(KML_Complex8), value :: a
        type(KML_Complex8), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_caxpyi    

    integer(C_INT) function kml_sparse_zaxpyi(nz, a, x, indx, y) bind(C, name='kml_sparse_zaxpyi')
        import
        integer(C_INT), value :: nz
        type(KML_Complex16), value :: a
        type(KML_Complex16), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zaxpyi    

    integer(C_INT) function kml_sparse_sroti(nz, x, indx, y, c, s) bind(C, name='kml_sparse_sroti')
        import
        integer(C_INT), value :: nz
        real(C_FLOAT), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_FLOAT), dimension(*), intent(out) :: y
        real(C_FLOAT), value :: c
        real(C_FLOAT), value :: s
    end function kml_sparse_sroti    

    integer(C_INT) function kml_sparse_droti(nz, x, indx, y, c, s) bind(C, name='kml_sparse_droti')
        import
        integer(C_INT), value :: nz
        real(C_DOUBLE), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_DOUBLE), dimension(*), intent(out) :: y
        real(C_DOUBLE), value :: c
        real(C_DOUBLE), value :: s
    end function kml_sparse_droti    

    integer(C_INT) function kml_sparse_sgthr(nz, y, x, indx) bind(C, name='kml_sparse_sgthr')
        import
        integer(C_INT), value :: nz
        real(C_FLOAT), dimension(*), intent(in) :: y
        real(C_FLOAT), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_sgthr    

    integer(C_INT) function kml_sparse_dgthr(nz, y, x, indx) bind(C, name='kml_sparse_dgthr')
        import
        integer(C_INT), value :: nz
        real(C_DOUBLE), dimension(*), intent(in) :: y
        real(C_DOUBLE), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_dgthr    

    integer(C_INT) function kml_sparse_cgthr(nz, y, x, indx) bind(C, name='kml_sparse_cgthr')
        import
        integer(C_INT), value :: nz
        type(KML_Complex8), dimension(*), intent(in) :: y
        type(KML_Complex8), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_cgthr    

    integer(C_INT) function kml_sparse_zgthr(nz, y, x, indx) bind(C, name='kml_sparse_zgthr')
        import
        integer(C_INT), value :: nz
        type(KML_Complex16), dimension(*), intent(in) :: y
        type(KML_Complex16), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_zgthr   

    integer(C_INT) function kml_sparse_sgthrz(nz, y, x, indx) bind(C, name='kml_sparse_sgthrz')
        import
        integer(C_INT), value :: nz
        real(C_FLOAT), dimension(*), intent(in) :: y
        real(C_FLOAT), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_sgthrz    

    integer(C_INT) function kml_sparse_dgthrz(nz, y, x, indx) bind(C, name='kml_sparse_dgthrz')
        import
        integer(C_INT), value :: nz
        real(C_DOUBLE), dimension(*), intent(in) :: y
        real(C_DOUBLE), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_dgthrz    

    integer(C_INT) function kml_sparse_cgthrz(nz, y, x, indx) bind(C, name='kml_sparse_cgthrz')
        import
        integer(C_INT), value :: nz
        type(KML_Complex8), dimension(*), intent(in) :: y
        type(KML_Complex8), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_cgthrz    

    integer(C_INT) function kml_sparse_zgthrz(nz, y, x, indx) bind(C, name='kml_sparse_zgthrz')
        import
        integer(C_INT), value :: nz
        type(KML_Complex16), dimension(*), intent(in) :: y
        type(KML_Complex16), dimension(*), intent(out) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
    end function kml_sparse_zgthrz   

    integer(C_INT) function kml_sparse_ssctr(nz, x, indx, y) bind(C, name='kml_sparse_ssctr')
        import
        integer(C_INT), value :: nz
        real(C_FLOAT), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_ssctr    

    integer(C_INT) function kml_sparse_dsctr(nz, x, indx, y) bind(C, name='kml_sparse_dsctr')
        import
        integer(C_INT), value :: nz
        real(C_DOUBLE), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dsctr    

    integer(C_INT) function kml_sparse_csctr(nz, x, indx, y) bind(C, name='kml_sparse_csctr')
        import
        integer(C_INT), value :: nz
        type(KML_Complex8), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_csctr    

    integer(C_INT) function kml_sparse_zsctr(nz, x, indx, y) bind(C, name='kml_sparse_zsctr')
        import
        integer(C_INT), value :: nz
        type(KML_Complex16), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zsctr   

    integer(C_INT) function kml_sparse_sdoti(nz, x, indx, y, doti) bind(C, name='kml_sparse_sdoti')
        import
        integer(C_INT), value :: nz
        real(C_FLOAT), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_FLOAT), dimension(*), intent(in) :: y
        type(C_PTR), value :: doti
    end function kml_sparse_sdoti    

    integer(C_INT) function kml_sparse_ddoti(nz, x, indx, y, doti) bind(C, name='kml_sparse_ddoti')
        import
        integer(C_INT), value :: nz
        real(C_DOUBLE), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        real(C_DOUBLE), dimension(*), intent(in) :: y
        type(C_PTR), value :: doti
    end function kml_sparse_ddoti    

    integer(C_INT) function kml_sparse_cdotci_sub(nz, x, indx, y, dotci) & 
            bind(C, name='kml_sparse_cdotci_sub')
        import
        integer(C_INT), value :: nz
        type(KML_Complex8), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex8), dimension(*), intent(in) :: y
        type(C_PTR), value :: dotci
    end function kml_sparse_cdotci_sub    

    integer(C_INT) function kml_sparse_zdotci_sub(nz, x, indx, y, dotci) & 
            bind(C, name='kml_sparse_zdotci_sub')
        import
        integer(C_INT), value :: nz
        type(KML_Complex16), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex16), dimension(*), intent(in) :: y
        type(C_PTR), value :: dotci
    end function kml_sparse_zdotci_sub 

    integer(C_INT) function kml_sparse_cdotui_sub(nz, x, indx, y, dotui) & 
            bind(C, name='kml_sparse_cdotui_sub')
        import
        integer(C_INT), value :: nz
        type(KML_Complex8), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex8), dimension(*), intent(in) :: y
        type(C_PTR), value :: dotui
    end function kml_sparse_cdotui_sub    

    integer(C_INT) function kml_sparse_zdotui_sub(nz, x, indx, y, dotui) & 
            bind(C, name='kml_sparse_zdotui_sub')
        import
        integer(C_INT), value :: nz
        type(KML_Complex16), dimension(*), intent(in) :: x
        integer(C_INT), dimension(*), intent(in) :: indx
        type(KML_Complex16), dimension(*), intent(in) :: y
        type(C_PTR), value :: dotui
    end function kml_sparse_zdotui_sub 

    integer(C_INT) function kml_set_thread_num(thread_num) bind(C, name='kml_set_thread_num')
        import
        integer(C_INT), value :: thread_num
    end function kml_set_thread_num

    integer(C_INT) function kml_sparse_scsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y) bind(C, name='kml_sparse_scsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scsrtrsv

    integer(C_INT) function kml_sparse_dcsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y) bind(C, name='kml_sparse_dcsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcsrtrsv

    integer(C_INT) function kml_sparse_ccsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y) bind(C, name='kml_sparse_ccsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccsrtrsv

    integer(C_INT) function kml_sparse_zcsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y) bind(C, name='kml_sparse_zcsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcsrtrsv
    
    integer(C_INT) function kml_csparse_scsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y) bind(C, name='kml_csparse_scsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        real(C_FLOAT), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_csparse_scsrtrsv

    integer(C_INT) function kml_csparse_dcsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y)bind(C, name='kml_csparse_dcsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        real(C_DOUBLE), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_csparse_dcsrtrsv

    integer(C_INT) function kml_csparse_ccsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y)bind(C, name='kml_csparse_ccsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        type(KML_Complex8), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_csparse_ccsrtrsv

    integer(C_INT) function kml_csparse_zcsrtrsv( &
            uplo, opt, diag, m, a, ia, ja, x, y)bind(C, name='kml_csparse_zcsrtrsv')
        import
        integer(C_INT), value :: uplo
        integer(C_INT), value :: opt
        integer(C_INT), value :: diag
        integer(C_INT), value :: m
        type(KML_Complex16), dimension(*), intent(in) :: a
        integer(C_INT), dimension(*), intent(in) :: ia
        integer(C_INT), dimension(*), intent(in) :: ja
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_csparse_zcsrtrsv

    integer(C_INT) function kml_sparse_scscmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) bind(C, name='kml_sparse_scscmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_FLOAT), value :: beta
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_scscmm

    integer(C_INT) function kml_sparse_dcscmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) bind(C, name='kml_sparse_dcscmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_DOUBLE), value :: beta
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_dcscmm

    integer(C_INT) function kml_sparse_ccscmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) bind(C, name='kml_sparse_ccscmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex8), value :: beta
        type(KML_Complex8), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_ccscmm

    integer(C_INT) function kml_sparse_zcscmm( &
            opt, m, n, k, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, beta, c, ldc) bind(C, name='kml_sparse_zcscmm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        integer(C_INT), value :: k
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex16), value :: beta
        type(KML_Complex16), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_zcscmm
    
    integer(C_INT) function kml_sparse_scscmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_scscmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), value :: beta
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scscmv
    
    integer(C_INT) function kml_sparse_dcscmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_dcscmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), value :: beta
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcscmv

    integer(C_INT) function kml_sparse_ccscmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_ccscmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), value :: beta
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccscmv

    integer(C_INT) function kml_sparse_zcscmv( &
            opt, m, k, alpha, matdescra, val, indx, pntrb, pntre, x, beta, y) bind(C, name='kml_sparse_zcscmv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: k
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), value :: beta
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcscmv

    integer(C_INT) function kml_sparse_scsrsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_scsrsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_scsrsm
    
    integer(C_INT) function kml_sparse_dcsrsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_dcsrsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_dcsrsm

    integer(C_INT) function kml_sparse_ccsrsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_ccsrsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex8), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_ccsrsm

    integer(C_INT) function kml_sparse_zcsrsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_zcsrsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex16), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_zcsrsm

    
    integer(C_INT) function kml_sparse_scscsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_scscsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_FLOAT), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_scscsm
    
    integer(C_INT) function kml_sparse_dcscsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_dcscsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        real(C_DOUBLE), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_dcscsm

    integer(C_INT) function kml_sparse_ccscsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_ccscsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex8), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_ccscsm

    integer(C_INT) function kml_sparse_zcscsm( &
            opt, m, n, alpha, matdescra, val, indx, pntrb, pntre, b, ldb, c, ldc) bind(C, name='kml_sparse_zcscsm')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        integer(C_INT), value :: n
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: b
        integer(C_INT), value :: ldb
        type(KML_Complex16), dimension(*), intent(out) :: c
        integer(C_INT), value :: ldc
    end function kml_sparse_zcscsm
    
    integer(C_INT) function kml_sparse_scscsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_scscsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_FLOAT), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_FLOAT), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_FLOAT), dimension(*), intent(in) :: x
        real(C_FLOAT), dimension(*), intent(out) :: y
    end function kml_sparse_scscsv

    integer(C_INT) function kml_sparse_dcscsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_dcscsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        real(C_DOUBLE), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        real(C_DOUBLE), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        real(C_DOUBLE), dimension(*), intent(in) :: x
        real(C_DOUBLE), dimension(*), intent(out) :: y
    end function kml_sparse_dcscsv

    integer(C_INT) function kml_sparse_ccscsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_ccscsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex8), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex8), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex8), dimension(*), intent(in) :: x
        type(KML_Complex8), dimension(*), intent(out) :: y
    end function kml_sparse_ccscsv

    integer(C_INT) function kml_sparse_zcscsv( &
            opt, m, alpha, matdescra, val, indx, pntrb, pntre, x, y) bind(C, name='kml_sparse_zcscsv')
        import
        integer(C_INT), value :: opt
        integer(C_INT), value :: m
        type(KML_Complex16), value :: alpha
        character(C_CHAR), dimension(*), intent(in) :: matdescra
        type(KML_Complex16), dimension(*), intent(in) :: val
        integer(C_INT), dimension(*), intent(in) :: indx
        integer(C_INT), dimension(*), intent(in) :: pntrb
        integer(C_INT), dimension(*), intent(in) :: pntre
        type(KML_Complex16), dimension(*), intent(in) :: x
        type(KML_Complex16), dimension(*), intent(out) :: y
    end function kml_sparse_zcscsv

end interface