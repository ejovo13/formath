program vec_test

use vector_m
use iso_fortran_env, only: real64, real32

implicit none

    integer, parameter :: dp = real64
    integer, parameter :: sp = real32

    type(vector) :: v1, v2, origin, v3, v4, u1, u2, u3, proj, v1_norm, v2_norm
    type(nd_vector(n = 10)) :: nd_10
    type(nd_vector(n = 3)) :: nd_3
    type(vector), dimension(2) :: orthonormal_basis
    type(vector), dimension(3) :: basis, copy_basis, ortho_3
    type(vector), dimension(2) :: basis_2
    real(dp) :: scalar

    v1 = vector(2)

    ! print *, "step 1 post vector_from"
    ! print *, "step 2 post vector_from"
    ! print *, "step 3 post vector_from"
    call v1%print_info()

    call v1%clear()

    v1 = [1, 2]

    call v1%print_info()
    ! print *, "v1.at(1) = ", v1%at(1)

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

    v2_norm = v2%normalize()
    v1_norm = v1%normalize()

    call v2_norm%print()
    call v1_norm%print()

    basis(1) = [1, 0, 0]
    basis(2) = [0, 1, 0]
    basis(3) = [0, 0, 1]

    copy_basis = basis

    call copy_basis(2)%set(2, 10)

    print *, "copy_basis: "
    call print_basis(copy_basis)




    print *, "original_basis: "
    call print_basis(basis)

    print *, "4 dimensional 0 basis of 2 vectors:"

    ! call basis_2%zero()

    ! call basis_2(1)%print()

    call print_basis(basis_2)
    ! call nd_test(nd_10)
    ! call nd_test(nd_3)

    orthonormal_basis = gram_schmidt(basis_2)
    ortho_3 = gram_schmidt(basis)

    call print_basis(orthonormal_basis)
    call print_basis(ortho_3)

    print *, "is orthonormal_basis orthonormal? ", is_orthonormal(orthonormal_basis)
    print *, "is ortho_3 orthonormal? ", is_orthonormal(ortho_3)


    

end program