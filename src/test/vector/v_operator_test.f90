program operator_test
!! Test the functional operators
use vector_m
use iso_fortran_env

implicit none

    type(vector) :: v1, v2, v3, res, z, x, y
    real(real64) :: scal_res

    v1 = [1, 2, 4]
    v2 = [6, 5, 10]
    z = [0, 0, 1]
    x = [1, 0, 0]
    y = [0, 1, 0]

    scal_res = v1 .dot. v2
    print *, "Dot product of ", v1%data(), " and ", v2%data(), " = ", scal_res

    v3 = v1 .proj. v2

    v2 = v2%normalized()
    v3 = v3%normalized()
    scal_res = v2 .dot. v3

    print*, "Dot product of v2 and v3 = ", scal_res

    res = v1 + v2

    v1 = [1, 5, 8]
    v2 = [1, -3, -5]
    res = v1 + v2
    print*, "v1: ", v1%data(), " + ", v2%data(), " = ", res%data()

    res = v1 - v2
    print*, "v1: ", v1%data(), " - ", v2%data(), " = ", res%data()

    res = v1 * 10
    print*, "v1: ", v1%data(), " * 10 = ", res%data()

    v2 = [100, 6, 5]
    res = v2 / 10
    print*, "v2: ", v2%data(), " /10 = ", res%data()

    res = v1 * v2
    print*, "v1: ", v1%data(), " * ", v2%data(), " = ", res%data()

    res = v1 / v2
    print*, "v1: ", v1%data(), " / ", v2%data(), " = ", res%data()

    res = 5 * v1
    print *, "5 * ", v1%data(), " = ", res%data()

    res = 5 / v1
    print *, "5 / ", v1%data(), " = ", res%data()

    res = -v1
    print *, "-(", v1%data(), ") = ", res%data()

    res = v1%orthogonalized(vector([0, 0, 1]))
    print *, "v1%orthogonalized([0, 0, 1]) = ", res%data()

    print *, "Scalar = " , (v1 .dot. z) / (z .dot. z)
    res = v1 .proj. z
    print *, "v1 .proj. z = ", res%data()

    res = v1%orthonormalized(vector([0, 0, 1]))
    print *, "v1%orthogonormalized([0, 0, 1]) = ", res%data()


    99 format(80("="))
    print 99
    print 99
    print *, "House holder transformation tests"
    print 99
    print 99

    res = v1 .hh. z
    print *, "v1 .hh. z = ", res%data()

    res = v1 .hh. y
    print *, "v1 .hh. y = ", res%data()

    res = v1 .hh. x
    print *, "v1 .hh. x = ", res%data()

    res = v1 .hhnorm. x
    ! res = res * -1/(res%at(1))
    print *, "v1 .hhnorm. x = ", res%data()

    res = v1 .hh. res
    print *, "v1 .hh. res = ", res%data()


    
end program
