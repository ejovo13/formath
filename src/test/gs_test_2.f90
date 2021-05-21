program gs_test_2

use iso_fortran_env
use matrix_m
use vector_m

implicit none

    type(matrix) :: m, ortho

    m = reshape([1, 2, 3, 4], [2, 2])

    print *, "M: "
    call m%print()

    ortho = m%gram_schmidt()

    print "(A)","Ortho: "
    call ortho%print()

end program