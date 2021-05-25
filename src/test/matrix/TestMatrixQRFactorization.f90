program mat_qr_test

use matrix_m
use vector_m
use iso_fortran_env, only: real64, real32

implicit none

    type(matrix) :: A, B, C, D, Q
    type(vector) :: norm
    integer :: ncol

    A = matrix([1, 2, 3, 4], 2, 2)

    call A%print()

    ncol = A%ncol()

    norm = A%get_col(1) .hhnorm. norm%id(dim=ncol, col=1)

    Q = Q%create_hh(norm)

    B = Q * A

    print *, "A = "
    call A%print()

    print *, "Times householder matrix: "
    call Q%print()

    print *, "Equals: "
    call B%print()


end program