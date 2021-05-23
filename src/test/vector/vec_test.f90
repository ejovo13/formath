program vec_test

use vector_m
use matrix_m
use iso_fortran_env, only: real64, real32

implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32

    type(vector) :: v1, v2, origin, v3, v4, u1, u2, u3, proj, v1_norm, v2_norm
    ! type(nd_vector(n = 10)) :: nd_10
    ! type(nd_vector(n = 3)) :: nd_3
    type(vector), dimension(3) :: basis
    type(vector), dimension(2) :: basis_2
    type(matrix) :: test_matrix
    type(matrix) :: t2_mat, ortho_m
    type(vector) :: vec_ptr
    real(dp) :: scalar
    real(dp), allocatable, dimension(:) :: test_a
    real(dp), allocatable, dimension(:,:) :: reg_array

    integer, dimension(4, 4) :: A

    A = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], [4, 4])

    u1 = vector(3)
    u2 = vector(4)
    u3 = vector(5)

    v1 = vector(2)

    ! print *, "step 1 post vector_from"
    ! print *, "step 2 post vector_from"
    ! print *, "step 3 post vector_from"
    v1 = [1, 2]
    call v1%print_info()

    call v1%clear()

    
    
    call v1%print_info()
    ! print *, "v1.at(1) = ", v1%at(1)
    
    v1 = vector(3)
    v1 = [1, 3, 4]
    call v1%print_info()

    
    v2 = v1
    print *, "Called v2 = v1"

    call v2%print_info()

    scalar = v2 .dot. v1

    print *, "Dot product of v1 and v2: ", scalar

    print *, "sqrt(scalar) = ", sqrt(scalar)

    print *, "length of v1 = ", v1%length()

    print *, "v1%pnorm(1) = ", v1%pnorm(1)
    print *, "v1%pnorm(100) = ", v1%pnorm(100)
    print *, "v1%pnorm(200) = ", v1%pnorm(200)

    print *, "v1%at(10) = ", v1%at(3)

    ! ! v1 = v1 .dot. v1

    call v2%set(3, 10)

    call v1%print()
    call v2%print()
    call origin%zero(10)

    call origin%print()

    v3 = [0, 1, 0]
    v4 = [1, 0, 0]

    print *, "v3 .dot. v4 = ", v3 .dot. v4

    print *, "are v3 and v4 orthogonal?", v3%is_ortho(v4)
    print *, "are v1 and v2 orthogonal?", v2%is_ortho(v1)

    1 format(80("="))

    print *
    print *
    print 1
    print *, "Starting gram-schmidt testing"
    print 1
    print *
    print *

    v1 = [3, 1]
    v2 = [2, 2]

    basis_2(1) = v1
    basis_2(2) = v2

    call v1%print()
    call v2%print()

    print *, "v2 projected onto v1 = "

    print *, "v1 .dot. v2 = ", v1 .dot. v2
    print *, "v1 .dot. v1 = ", v1 .dot. v1

    proj = v2 .proj. v1
    call proj%print()

    print *, v2 .dot. v1

    v2_norm = v2%normalized()
    v1_norm = v1%normalized()

    call v2_norm%print()
    call v1_norm%print()

    basis(1) = [1, 0, 0]
    basis(2) = [0, 1, 0]
    basis(3) = [0, 0, 1]





    print *
    print *
    print 1
    print *, "Starting matrix testing"
    print 1
    print *
    print *

    call test_matrix%new(2, 3)

    vec_ptr = test_matrix%vec(1)

    call vec_ptr%print()

    call test_matrix%print()

    print *, "element(1, 1) = ", test_matrix%at(1, 1)


    t2_mat = A    
    test_matrix = t2_mat

    call test_matrix%set(1, 1, 10._real64)
    call test_matrix%set(3, 2, 10)

    
    print *, "test matrix: "
    call test_matrix%print()

    print *, "unmodified: "
    call t2_mat%print()

    ortho_m = test_matrix%gram_schmidt()

    print *, "Gram-schmidt solve:"

    call ortho_m%print()

    print *, "Is ortho_m orthonormal? ", ortho_m%is_orthonormal()

    test_matrix = reshape([1, 3, 4, 5, 9, -2], [2, 3])

    call test_matrix%print()

    test_a = v2%data()

    ortho_m = test_matrix%gram_schmidt()
    print *, "Gram-schmidt solution: "

    call ortho_m%print()

    reg_array = ortho_m%as_array()

    print *, reg_array

end program