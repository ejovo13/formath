program mat_ctor_test

use matrix_m
use vector_m

use iso_fortran_env, only: real64, real32

implicit none

    type(matrix) :: m1, m2, m3, m4, m5, m6

    m1 = matrix(reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [4,3]))
    m2 = matrix(real(reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [4,3]), real32))
    m3 = matrix(real(reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [4,3]), real64))

    m4 = matrix(3, 3, 13)

    m5 = matrix(m4)

    call m5%set_col(3, [1, 2, 3])

    m6 = matrix(3,3)


    call m1%print()
    call m2%print()
    call m3%print()
    call m4%print()
    call m5%print()
    call m6%print()


end program