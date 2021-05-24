program mat_qr_test

use matrix_m
use vector_m
use iso_fortran_env, only: real64, real32

implicit none

    type(matrix) :: A, B, C, D 

    A = matrix([1, 2, 3, 4], 4, 1)

    call A%print()


end program