program v_operator_sub_test

use vector_m

implicit none

    type(vector) :: v1, v2, v3

    v1 = [1, 2, 4]
    v2 = [6, 5, 10]
    v3 = v1

    call v3%proj(v1)
    print*, v3%data()

    print *, "Starting v1: ", v1%data()
    print *, "Starting v2: ", v2%data()

    call v1%plus(v2)
    print *, "v1%plus(v2) = ", v1%data()

    call v1%minus(v2)
    print *, "v1%minus(v2) = ", v1%data()

    call v1%times(v2)
    print *, "v1%times(v2) = ", v1%data()

    call v1%times(5)
    print *, "v1%times(5) = ", v1%data()

    call v1%div(100) 
    print *, "v1%div(100) = ", v1%data()

    call v1%div(v2)
    print *, "v1%div(v2) = ", v1%data()

    call v1%normalize()
    print *, "v1%normalize() = ", v1%data()

end program