program operator_test

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

    call v2%normalize()
    call v3%normalize()

    res = v1 + v2

end program
