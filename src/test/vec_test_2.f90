program vec_test_2

use matrix_m
use vector_m
use iso_fortran_env

implicit none

    type(vector) :: v1

    v1 = vector(2)

    v1 = [1, 2]

    call v1%print()

    v1 = [2]

    call v1%print_info()

    call test_sub(v1, [2._real64])

contains 


    subroutine test_sub(self, array) 

        class(vector) :: self
        real(real64), dimension(:) :: array

        if (size(array) /= self%n()) then
            error stop "Please don't compile me"
        end if

        print *, "Array with dim =", size(array), " passed"

    end subroutine

end program