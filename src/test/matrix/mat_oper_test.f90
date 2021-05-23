program mat_subop_test
    ! Test the function operators for matrices

use matrix_m
use vector_m
use iso_fortran_env, only: real64

implicit none

    type(matrix) :: m1, m2, m3, m4, m5, m6
    type(vector) :: v1, v2

    m1 = reshape([2,3,4,5], [2,2])
    m2 = reshape([1, 0, 0, 1], [2,2])

    v1 = [2, 4]

    m3 = m1 * m2

    call m1%print()
    print *, " times "
    call m2%print()
    print *, " equals "
    call m3%print()

    v2 = m1 * v1

    call m1%print()
    print *, "times"
    call v1%print()
    print *, " equals "
    call v2%print()
    
    m4 = m1 .o. m2

    call m4%print()
    print *

    m5 = m1 + m2

    call m1%print() 
    print *, "plus"
    call m2%print()
    print *, "equals"
    call m5%print()

    m6 = m4 - m3


    call m4%print() 
    print *, "minus"
    call m3%print()
    print *, "equals"
    call m6%print()




end program