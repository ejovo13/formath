program v_constructor_test
!! Test whether the assignments work as expected

use vector_m

implicit none

    type(vector) :: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13

    v1 = vector(5)                              !! Create empty vector(zeros) with 5 elements
    v2 = vector([1, 2, 3])                      !! Create dim 3 vector from [1, 2, 3]
    v3 = vector(4, 0)                           !! Create 4 dim vector and fill it with 0
    v4 = vector([1.0, 3.0, 4.0, 5.0, 3.0])
    v5 = vector(5, 1.6)    
    v6 = vector(10, 1.7d0)
    v7 = vector([1.4d0, -14.d0])
    v8 = vector(v7)

    v9 = v8%eye()
    v10 = v10%eye(20)
    v11 = v11%eye(5, 3)
    v12 = v7%eye(col=1)
    v13 = v4%eye(dim=10, col=2)

    call v1%print_info()
    call v2%print_info()
    call v3%print_info()
    call v4%print_info()
    call v5%print_info()
    call v6%print_info()
    call v7%print_info()
    call v8%print()    
    call v9%print()    
    call v10%print()    
    call v11%print()    
    call v12%print()    
    call v13%print()    
    
end program