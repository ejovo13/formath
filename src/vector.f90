module vector_m
!! Attempt to play around with a vector object 

!! A vector has the following operations:

!!
!! Addition with another vector
!! Multiplication with a scalar 
!! Inner product
!! Norm (which norm?)
!! Projection of a onto b 
!!
!!

!! The internal representaion of a vector is an array of fixed size, pointed to by a pointer.
!! This way we can use the fundamental data of an array within our specific class
use iso_fortran_env, only: real64, real32, int32, int64

implicit none
private
public :: print_basis, gram_schmidt, is_orthonormal, deallocate_vector_data

real(real64), parameter :: vector_epsilon = 1d-14

type, public :: vector

    private
    real(real64), dimension(:), allocatable :: v
    integer :: dim = 0

contains 

    private
    procedure, public :: size => vector_size
    procedure, public :: new => new_constructor
    procedure, public :: clear => clear_vector
    !! Deallocate the data, set v_allocated to false, set dim to 0
    procedure, public :: print_info => vector_print_info
    procedure, public :: print => vector_print_coordinates
    procedure, public :: at => vector_at_index
    generic, public :: set => set_int_, set_r32_, set_r64_
    procedure, public :: length => vector_euclidiean_norm
    procedure, public :: pnorm => vector_pnorm
    procedure, public :: normalize => vector_normalize
    procedure, public :: is_ortho => vector_is_orthogonal
    procedure, public :: is_normal => vector_is_normal
    procedure, public :: as_array => vector_as_array

    procedure, public :: zero => vector_zero
    
    procedure, public :: allocated => vector_is_allocated


    procedure, public :: alloc_ => allocate_vector_data
    procedure, public :: dealloc_ => deallocate_vector_data
    procedure :: dot_ => vector_dot_vector
    procedure :: proj_ => vector_proj_vector
    procedure :: conform_ => vector_conform
    procedure :: scalar_mult_int_ => vector_times_scalar_int
    procedure :: scalar_mult_r32_ => vector_times_scalar_r32
    procedure :: scalar_mult_r64_ => vector_times_scalar_r64
    procedure :: scalar_div_int_ => vector_div_scalar_int
    procedure :: scalar_div_r32_ => vector_div_scalar_r32
    procedure :: scalar_div_r64_ => vector_div_scalar_r64
    procedure :: set_int_ => vector_set_index_int
    procedure :: set_r32_ => vector_set_index_r32
    procedure :: set_r64_ => vector_set_index_r64
    procedure :: minus_ => vector_minus_vector
    procedure :: plus_ => vector_plus_vector


    generic, public :: assignment(=) => from_array_int_, from_array_r32_, from_array_r64_, from_vector_, from_int_, &
                                        from_r32_, from_r64_
    generic, public :: operator(.dot.) => dot_
    generic, public :: operator(.proj.) => proj_
    generic, public :: operator(*) => scalar_mult_int_, scalar_mult_r32_, scalar_mult_r64_
    generic, public :: operator(/) => scalar_div_int_, scalar_div_r32_, scalar_div_r64_
    generic, public :: operator(+) => plus_
    generic, public :: operator(-) => minus_


    procedure :: from_int_ => vector_from_int
    procedure :: from_r32_ => vector_from_r32
    procedure :: from_r64_ => vector_from_r64
    procedure :: from_array_int_ => vector_from_array_int
    procedure :: from_array_r32_ => vector_from_array_r32
    procedure :: from_array_r64_ => vector_from_array_r64
    procedure :: from_vector_ => vector_from_vector

    final :: vector_destructor    

end type
interface vector
    
    procedure :: vector_constructor_int
    procedure :: vector_constructor_r32
    procedure :: vector_constructor_r64
    procedure :: vector_constructor_dim
    procedure :: vector_constructor_dim_value_int
    procedure :: vector_constructor_dim_value_r32
    procedure :: vector_constructor_dim_value_r64


end interface 

contains 


!=============================================================================!
!=                               Constructors                                =!
!=============================================================================!

    pure subroutine new_constructor(self, dim) 
    !! allocate the proper space for our vector and set the dimension
        class(vector), intent(inout) :: self
        integer, intent(in) :: dim

        self%dim = dim
        allocate(self%v(dim))
        

    end subroutine

    pure function vector_constructor_int(array) result(this)
    
        integer, dimension(:), intent(in) :: array
        type(vector) :: this

        call this%new(size(array))
        this%v = array

    end function

    pure function vector_constructor_r32(array) result(this)
    
        real(real32), dimension(:), intent(in) :: array
        type(vector) :: this

        call this%new(size(array))
        this%v = array

    end function

    pure function vector_constructor_r64(array) result(this)
    
        real(real64), dimension(:), intent(in) :: array
        type(vector) :: this

        call this%new(size(array))
        this%v = array

    end function

    elemental function vector_constructor_dim(dim) result(this)

        integer, intent(in) :: dim
        type(vector) :: this

        call this%new(dim)
        this%v = 0

    end function

    elemental function vector_constructor_dim_value_int(dim, val) result(this)

        integer, intent(in) :: dim
        integer, intent(in) :: val

        type(vector) :: this

        call this%new(dim)
        this%v = val

    end function

    elemental function vector_constructor_dim_value_r32(dim, val) result(this)

        integer, intent(in) :: dim
        real(real32), intent(in) :: val

        type(vector) :: this

        call this%new(dim)
        this%v = val

    end function

    elemental function vector_constructor_dim_value_r64(dim, val) result(this)

        integer, intent(in) :: dim
        real(real64), intent(in) :: val

        type(vector) :: this

        call this%new(dim)
        this%v = val

    end function

    elemental subroutine clear_vector(self)

        class(vector), intent(inout) :: self

        if (self%allocated()) then
            call self%dealloc_()
        else
            self%dim = 0
        end if

    end subroutine

    elemental subroutine allocate_vector_data(self, dim) 

        class(vector), intent(inout) :: self
        integer, intent(in) :: dim

        integer :: ierr

        allocate(self%v(dim), STAT=ierr)
        self%dim = dim

        if (ierr /= 0) error stop "Error allocating vector"

    end subroutine

    elemental subroutine deallocate_vector_data(self) 

        class(vector), intent(inout) :: self
        integer :: ierr

        deallocate(self%v, STAT=ierr)
        self%dim = 0

        if (ierr /= 0) error stop "Error allocating vector"

    end subroutine

!=============================================================================!
!=                         Assigment Functions                               =!
!=============================================================================!
    pure subroutine vector_from_int(self, array) 

        class(vector), intent(inout) :: self
        integer, intent(in) :: array

        self%v = array ! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_r32(self, array) 

        class(vector), intent(inout) :: self
        real(real32), intent(in) :: array

        self%v = array! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_r64(self, array) 

        class(vector), intent(inout) :: self
        real(real64), intent(in) :: array

        self%v = array ! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_array_int(self, array) 

        class(vector), intent(inout) :: self
        integer, dimension(:), intent(in) :: array

        self%dim = size(array)
        self%v = array ! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_array_r32(self, array) 

        class(vector), intent(inout) :: self
        real(real32), dimension(:), intent(in) :: array

        self%dim = size(array)
        self%v = array! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_array_r64(self, array) 

        class(vector), intent(inout) :: self
        real(real64), dimension(:), intent(in) :: array

        self%dim = size(array)
        self%v = array ! Copy the contents of array into self

    end subroutine

    elemental subroutine vector_from_vector(self, v1)

        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v1

        self%dim = v1%dim
        self%v = v1%v ! Copy the array contents

    end subroutine

!=============================================================================!
!=                            State Functions                                =!
!=============================================================================!

    elemental function vector_conform(self, v2) result(bool)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        logical :: bool

        bool = (self%dim == v2%dim)

    end function

    elemental function vector_is_orthogonal(self, v2, eps) result(bool)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2
        real(real64), intent(in), optional :: eps

        logical :: bool

        if(present(eps)) then

            bool = abs(self .dot. v2) < eps

        else

            bool = abs(self .dot. v2) < vector_epsilon

        end if

    end function

    elemental function vector_is_normal(self) result(bool)

        class(vector), intent(in) :: self
        logical :: bool

        if (self%length() - 1 < vector_epsilon) then
            bool = .true.
        else 
            bool = .false.
        end if

    end function

    elemental function vector_is_allocated(self) result(bool)

        class(vector), intent(in) :: self

        logical :: bool

        if (allocated(self%v)) then
            bool = .true.
        else
            bool = .false.
        end if

    end function

    elemental function vector_size(self) result(n)

        class(vector), intent(in) :: self
        integer :: n

        n = self%dim

    end function

    
!=============================================================================!
!=                            Print Functions                                =!
!=============================================================================!
    
    subroutine vector_print_info(self) 
        
        class(vector), intent(in) :: self
        
        if (.not. self%allocated()) then
            write(*,*) "Vector not allocated"
        else 
            write(*,*) "dimension: ", self%dim
            write(*,*) "data: ", self%v
            write(*,*) "allocated: ", self%allocated()
        end if
        
    end subroutine
    
    subroutine vector_print_coordinates(self) 
        
        class(vector), intent(in) :: self
        
        print *, self%v
        
    end subroutine
    
!=============================================================================!
!=                           Access Functions                                =!
!=============================================================================!
    
    elemental function vector_at_index(self, index) result(x_n)
        
        class(vector), intent(in) :: self
        integer, intent(in) :: index
        
        real(real64) :: x_n
        
        if (index > self%dim) error stop "out of bounds error"
        
        x_n = self%v(index)
        
    end function

    elemental subroutine vector_set_index_int(self, index, val)
        
        class(vector), intent(inout) :: self
        integer, intent(in) :: index
        integer, intent(in) :: val
        
        if (index > self%dim) error stop "out of bounds error"
        
        self%v(index) = real(val, real64)
        
    end subroutine

    elemental subroutine vector_set_index_r32(self, index, val)
        
        class(vector), intent(inout) :: self
        integer, intent(in) :: index
        real(real32), intent(in) :: val
        
        if (index > self%dim) error stop "out of bounds error"
        
        self%v(index) = real(val, real64)
        
    end subroutine

    elemental subroutine vector_set_index_r64(self, index, val)
        
        class(vector), intent(inout) :: self
        integer, intent(in) :: index
        real(real64), intent(in) :: val
        
        if (index > self%dim) error stop "out of bounds error"
        
        self%v(index) = val
        
    end subroutine

    pure function vector_as_array(self) result(array)

        class(vector), intent(in) :: self
        real(real64), dimension(self%dim) :: array

        array = self%v

    end function
    
!=============================================================================!
!=                             Norm Functions                                =!
!=============================================================================!
    
    elemental function vector_euclidiean_norm(self) result(length)
        
        class(vector), intent(in) :: self
        real(real64) :: length
        
        length = self%pnorm(2)
        
    end function
    
    elemental function vector_pnorm(self, p) result(pnorm)
        
        class(vector), intent(in) :: self
        integer, intent(in) :: p
        real(real64) :: power, pnorm
        
        if (p < 1) error stop "P-norm must be greater than or equal to 1"
        
        power = 1._real64 / p
        
        pnorm = sum(self%v**p)**power
        
    end function
    
    elemental function vector_normalize(self) result(normalized_vector)
    !! Normalize a vector such that its euclidian norm is 1

        class(vector), intent(in) :: self        
        type(vector) :: normalized_vector

        real(real64) :: norm

        norm = self%length()

        normalized_vector = (self/norm)

    end function

!=============================================================================!
!=                            Operator Functions                             =!
!=============================================================================!
    
    elemental function vector_dot_vector(self, v2) result(inner_product)
        !! Calculate the inner product of two vectors
        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2
        
        real(real64) :: inner_product
        
        if (self%conform_(v2)) then
            
            inner_product = sum(self%v * v2%v)
            
        else 
            
            error stop "Cannot take the inner product of nonconforming vectors"

        end if
        
    end function
    
    elemental function vector_proj_vector(self, v2) result(v3)
    !! Project vector self ONTO v2 \(proj_v2(self)\)
        
        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2
        type(vector) :: v3

        real(real64) :: scalar

        scalar = (self .dot. v2) / (v2 .dot. v2)
        
        v3 = v2 * scalar
        ! Allocate the space for a new vector
        
    end function

    elemental function vector_plus_vector(self, v2) result(v3)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        if (self%conform_(v2)) then
            call v3%new(self%dim)
            v3%v = self%v! + v2%v
        else 
            error stop "Cannot add nonconforming vectors"
        end if

    end function

    elemental function vector_minus_vector(self, v2) result(v3)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        if (self%conform_(v2)) then
            call v3%new(self%dim)
            v3%v = self%v - v2%v
        else 
            error stop "Cannot add nonconforming vectors"
        end if

    end function

    elemental function vector_times_scalar_int(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        integer, intent(in) :: scalar 
        type(vector) :: v2

        call v2%new(self%dim)

        v2%v = self%v * scalar

    end function

    elemental function vector_times_scalar_r32(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real32), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new(self%dim)

        v2%v = self%v * scalar
    
    end function

    elemental function vector_times_scalar_r64(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real64), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new(self%dim)

        v2%v = self%v * scalar
    
    end function

    elemental function vector_div_scalar_int(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        integer, intent(in) :: scalar 
        type(vector) :: v2

        call v2%new(self%dim)

        v2%v = self%v / scalar

    end function

    elemental function vector_div_scalar_r32(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real32), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new(self%dim)

        v2%v = self%v / scalar
    
    end function

    elemental function vector_div_scalar_r64(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real64), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new(self%dim)

        v2%v = self%v / scalar
    
    end function

!=============================================================================!
!=                           Static Functions                                =!
!=============================================================================!

    elemental subroutine vector_zero(self, dim)

        class(vector), intent(inout) :: self
        integer, intent(in), optional :: dim

        if (present(dim)) then
            call self%new(dim)
        else
            call self%new(1)
        end if

        self%v = 0

    end subroutine



!=============================================================================!
!=                                Destructor                                 =!
!=============================================================================!    
    
    subroutine vector_destructor(self) 
    
        type(vector), intent(inout) :: self
    
        write(*,*) "vector_destructor called"
    
        if(self%allocated()) then
            print *, "deallocating fields"
            deallocate(self%v)
        end if
    
    end subroutine
    
!=============================================================================!
!=                              Gram-Schmidt                                 =!
!=============================================================================!
    
    function gram_schmidt(vector_array) result(ortho)

        class(vector), dimension(:) :: vector_array
        type(vector), dimension(:), allocatable :: ortho

        integer i, j, k, n

        n = vector_array(1)%dim
        k = size(vector_array)

        allocate(ortho(k))
        call ortho%zero(n)

        print *, "dimension of input basis = ", vector_array(1)%dim
        print *, "number of basis vectors = ", size(vector_array)
        print *, "orthonormal_basis set to 0"

        ortho(1) = vector_array(1)%normalize()

        do i = 2,k

            ortho(i) = vector_array(i)

            do j = 1,i-1

                ortho(i) = ortho(i) - (ortho(i) .proj. ortho(j))

            end do

            ortho(i) = ortho(i)%normalize()

        end do

        
    end function

    subroutine print_basis(basis) 

        class(vector), dimension(:), intent(in) :: basis

        integer :: i, k
        
        k = size(basis)

        do i = 1, k

            call basis(i)%print()

        end do

    end subroutine 
    
    function is_orthonormal(basis) result(bool)

        class(vector), dimension(:), intent(in) :: basis
        logical :: bool

        integer :: k, i, j

        k = size(basis)

        if(all(basis%is_normal())) then

            do i = 1,k
                do j = i+1,k
                    if (.not. basis(i)%is_ortho(basis(j))) then
                        bool = .false.
                        return
                    end if
                end do
            end do

            bool = .true.

        else 

            bool = .false.

        end if     

    end function
    
end module