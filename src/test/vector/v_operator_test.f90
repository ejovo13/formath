program operator_test
!! Test the functional operators
use vector_m
use iso_fortran_env

implicit none

    type(vector) :: v1, v2, v3, res
    real(real64) :: scal_res

    v1 = [1, 2, 4]
    v2 = [6, 5, 10]
    
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








end program
