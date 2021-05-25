program v_assignment_test
!! Test whether the assignments work as expected

use vector_m

implicit none

    type(vector) :: v1, v2, v3, v4, v5, v6, v7

    v1 = 4
    v2 = [4, 2]
    v3 = [1.0, 2.0]
    v4 = [3.d0, 5.7d0, 109.d0, -50.d0]
    v5 = v4
    v6 = v5
    v6 = 15.0
    v7 = vector(20)
    v7 = -14.d0

    call v1%print_info()
    call v2%print_info()
    call v3%print_info()
    call v4%print_info()
    call v5%print_info()
    call v6%print_info()
    call v7%print_info()

end program