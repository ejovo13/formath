program v_n_constructor_test
    !! Test whether assigment or construction is faster
    use vector_m
    
    implicit none
    
        type(vector) :: v1
    
        integer :: i
        integer, parameter :: N_LOOPS = 1E7
    
        do i = 1, N_LOOPS
    
            v1 = vector([1, 3, 4, 5])
    
        end do
    
        call v1%print_info()
    
    end program