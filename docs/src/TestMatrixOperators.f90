program mat_subop_test
    ! Test the function operators for matrices

use matrix_m
use vector_m
use iso_fortran_env, only: real64

implicit none

    type(matrix) :: m1, m2, m3, m4, m5, m6, i1, m7, hh, m8, m9, o1, o2, m10, m11
    type(vector) :: v1, v2, v3, v4, v5, v_temp
    real(real64), dimension(:,:), allocatable :: temp_a

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

    i1 = i1%id(10)

    call i1%print()

    print *, "V1 = "
    call v1%print()
    print *, "V2 = "
    call v2%print()

    m7 = v1 .outer. v2

    print *, "Outer product: "
    call m7%print()

    v3 = [1, 2, 3]
    v4 = [1, 0, 0]

    v5 = v3 .hh. v4

    print *, v3%data(), " rotated over ", v4%data(), " = ", v5%data()

       
    hh = hh%create_hh(v4)
    call hh%print()

    hh = hh%create_hh(vector([0, 1, 0]))
    call hh%print()
    
    hh = hh%create_hh(vector([0, 0, 1]))
    call hh%print()

    
    v5 = hh * v3
    
    call v5%print()
    
    hh = hh%create_hh(vector([0, 0, 0, 0, 1]))
    call hh%print()

    print *
    print *
    print *, "Matrix decomposition methods"
    print *
    print *

    allocate(temp_a(3,2))

    call random_number(temp_a)

    m8 = temp_a

    call m8%print()

    o1 = m8%gram_schmidt()

    print *, "Calling gram_schmidt"

    print *
    print *

    call o1%print()

    print *, "is o1 orthonormal? : ", o1%is_orthonormal()


    v_temp = o1%get_col(1)

    print *, "Col 1 of o1 length = ", v_temp%length()

    deallocate(temp_a)
    allocate(temp_a(2, 4))

    call random_number(temp_a)

    m9 = temp_a

    o2 = m9%gram_schmidt()

    print *
    call m9%print()

    print*

    call o2%print()
    print *

    m10 = reshape([1, 2, 2, 4], [2,2])

    m11 = m10 ** 0

    call m10%print()

    print*, "m10 ** 0 = "
    call m11%print()

    m11 = m10**1

    print *, "m10 ** 1 = "
    call m11%print()

    m11 = m10**2

    print *, "m10 ** 2 = "
    call m11%print()

    m11 = m10**3

    print *, "m10 ** 3 = "
    call m11%print()

    m11 = m10**4

    print *, "m10 ** 4 = "
    call m11%print()


    
    
end program