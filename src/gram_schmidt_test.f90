program gs_test

use vector_m
use matrix_m
use iso_fortran_env, only: real64, real32

implicit none

    integer :: num_loops = 0, i
    character(20) :: num_loops_command_line

    real(real64), dimension(:), allocatable :: x
    type(matrix) :: base, ortho

    call get_command_argument(1, num_loops_command_line)

    read(num_loops_command_line, "(I10)") num_loops

    base = reshape([61, 27, 66, 69, 75, 46, 9, 23, 92, 16, 83, 54, 100, 8, 45, 11, 97, 1, 78, 82, 87, 9, 40, 26, 81], [5, 5])

    call base%print()

    do i = 1, num_loops

        ortho = base%gram_schmidt()

    end do

    call ortho%print()

    x = [1, 2]

    print *, x

end program