program v_op_speed

use vector_m

implicit none

    integer, parameter :: N_LOOPS = 1E8
    integer :: i
    type(vector) :: v1

    v1 = [1, 5, -4, 10]

    do i = 1, N_LOOPS

        v1 = v1 * 5

    end do

end program