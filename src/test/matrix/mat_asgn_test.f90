program m_asgn_test

use matrix_m
use vector_m
use iso_fortran_env, only: real64

implicit none

    type(matrix) :: m1, m2, m3, m4
    type(vector) :: v1

    m1 = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9],[3, 2])
    m2 = m1
    m3 = m2%as_array()
    call m3%set_row(2, [10, 20])
    call m3%set_col(2, [100, 200, 300])

    v1 = m3%get_col(1)

    call v1%print()

    v1 = m3%get_col(2)

    call v1%print()


    m4 = real(m3%as_array(), real64)
    call m4%set_row(3, [50, 60])

    call m1%print()
    call m2%print()
    call m3%print()
    call m4%print()



end program