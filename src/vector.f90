module vector_m
!! Attempt to play around with a vector object 

!! A vector, \(v\) has the following operations:

!!
!! Addition with another vector
!! Multiplication with a scalar 
!! Inner product
!! Norm (which norm?)
!! Projection of a onto b 
!!
!!

! // TODO Add support for array operations

!! The internal representaion of a vector is an allocatable array.
!! This way we can use the fundamental data of an array within our specific class
use iso_fortran_env, only: real64, real32, int32, int64

implicit none
private
public :: print_basis, gram_schmidt, is_orthonormal, operator(*), operator(/)!, deallocate_vector_data

real(real64), parameter :: vector_epsilon = 1d-14

type, public :: vector

    private
    real(real64), dimension(:), allocatable :: v
    integer :: dim = 0

contains 

    private
    procedure, public :: size => vector_size
    !! Return the number of elements of the currently allocated vector
    procedure, public :: clear => clear_vector
    !! Deallocate the data, set dim to 0
    procedure, public :: print_info => vector_print_info
    !! Print diagnostic information about the vector
    procedure, public :: print => vector_print_coordinates
    !! Print only the data stored in the vector object as a row vector
    procedure, public :: at => vector_at_index
    !! Return the element \(x_i\) 
    generic, public :: set => set_int_, set_r32_, set_r64_
    !! Set the value of the element \(x_i\) at index \(i\)
    procedure, public :: length => vector_euclidiean_norm
    !! Calculate the euclidean norm of a vector \(v\)
    procedure, public :: pnorm => vector_pnorm
    !! Calculate the pnorm of a vector \(v\)
    procedure, public :: normalize => vector_normalize
    !! Normalize the elements of the passed vector \(v\)
    !! Normaliz**e** is a **subroutine** such that it alters the elements of the passed vector \(v\) to avoid costs of copying involved with a function
    procedure, public :: normalized => vector_normalized
    !! Return a normalized vector \(n\) pointing in the same direction as \(v\)
    !!@Note
    !! The function normaliz**ed** is a **function** such that it returns a normalized version of the passed vector \(v\)
    procedure, public :: orthogonalize => vector_orthogonalize
    !! Orthogonalize a vector \(v\) against a passed **normalized** vector \(n\)
    !!@Note
    !!A future version may just check if the passed vector is normalized by testing a "normalized" logical type that will be stored in a vector.
    procedure, public :: orthogonalized => vector_orthogonalized
    !! Return a vector \(v\) that is orthogonalized against a passed **normalized** vector \(n\)
    !!@Note
    !! This is a function that returns a new vector \(v\)
    procedure, public :: orthonormalize => vector_orthonormalize
    !! Orthogonalize and normalize a vector \(v\) against a passed **normalized** vector \(n\)
    !!@Note
    !! This is a subroutine that modifies the passed vector \(v\)
    procedure, public :: orthonormalized => vector_orthonormalized
    !! Return an orthogonalized and normalized vector \(v\) against a passed **normalized** vector \(n\)
    !!@Note
    !! This is a function that returns a new vector \(v\)
    procedure, public :: householder_transform => vector_householder_sub
    !! Rotate a passed vector \(v\) about the hyper plane described by the passed **normalized** vector \(n\)
    !!@Note
    !! This is a subroutine that modifies the passed vector  
    
    procedure, public :: eye => vector_constructor_eye
    procedure, public :: is_ortho => vector_is_orthogonal
    procedure, public :: is_normal => vector_is_normal
    procedure, public :: data => vector_as_array
    
    procedure, public :: zero => vector_zero
    
    procedure, public :: allocated => vector_is_allocated
    
    
    procedure, public :: alloc_ => allocate_vector_data
    procedure, public :: dealloc_ => deallocate_vector_data
    procedure :: new_ => new_constructor
    procedure :: dot_ => vector_dot_vector
    procedure :: proj_ => vector_proj_vector
    procedure :: outer_ => vector_outer_vector
    procedure :: householder_ => vector_householder
    procedure :: householder_normal_ => vector_find_householder_normal
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
    procedure :: unary_minus_ => vector_unary_minus
    procedure :: plus_ => vector_plus_vector

    !===============Operators Functions==================!
    generic, public :: assignment(=) => from_array_int_, from_array_r32_, from_array_r64_, from_vector_, from_int_, &
                                        from_r32_, from_r64_
    generic, public :: operator(.dot.) => dot_
    generic, public :: operator(.inner.) => dot_
    generic, public :: operator(.outer.) => outer_
    generic, public :: operator(.proj.) => proj_
    generic, public :: operator(.hh.) => householder_
    generic, public :: operator(.hhnorm.) => householder_normal_
    generic, public :: operator(*) => scalar_mult_int_, scalar_mult_r32_, scalar_mult_r64_
    generic, public :: operator(.o.) => hadamard_vec_
    generic, public :: operator(/) => scalar_div_int_, scalar_div_r32_, scalar_div_r64_, div_vec_
    generic, public :: operator(+) => plus_
    generic, public :: operator(-) => minus_, unary_minus_

    !===============Operator Subroutines==================!
    ! The point of an operator subroutine is to alter the passed object. This cuts down on copying
    ! the data between functions
    generic, public :: times => times_int_sub_, times_r32_sub_, times_r64_sub_, times_vec_sub_
    generic, public :: div => div_int_sub_, div_r32_sub_, div_r64_sub_, div_vec_sub_
    generic, public :: proj => project_onto_sub_
    generic, public :: plus => plus_vector_sub_
    generic, public :: minus => minus_vector_sub_

    procedure :: times_int_sub_ => vector_times_scalar_int_sub
    procedure :: times_r32_sub_ => vector_times_scalar_r32_sub
    procedure :: times_r64_sub_ => vector_times_scalar_r64_sub
    procedure :: times_vec_sub_ => vector_times_vector_sub
    procedure :: div_int_sub_ => vector_div_scalar_int_sub
    procedure :: div_r32_sub_ => vector_div_scalar_r32_sub
    procedure :: div_r64_sub_ => vector_div_scalar_r64_sub
    procedure :: div_vec_sub_ => vector_div_vector_sub
    procedure :: project_onto_sub_ => vector_proj_vector_sub
    procedure :: plus_vector_sub_ => vector_plus_vector_sub
    procedure :: minus_vector_sub_ => vector_minus_vector_sub



    procedure :: from_int_ => vector_from_int
    procedure :: from_r32_ => vector_from_r32
    procedure :: from_r64_ => vector_from_r64
    procedure :: from_array_int_ => vector_from_array_int
    procedure :: from_array_r32_ => vector_from_array_r32
    procedure :: from_array_r64_ => vector_from_array_r64
    procedure :: from_vector_ => vector_from_vector
    procedure :: hadamard_vec_ => vector_hadamard_vector
    procedure :: div_vec_ => vector_div_vector

    final :: vector_destructor    

end type
interface vector
!! Construct a vector object <br>
!!@Note
!!A vector can be instantiated from integer and real data types but the underlying data will always
!!be stored using double precision.
!!
!! The following code snippet shows all of the valid ways to construct a new vector
!!```fortran
!!
!!type(vector) :: v1, v2, v3, v4, v5, v6, v7, v8
!!
!!v1 = vector([1, 2, 3])
!!v2 = vector([1.0, 2.0, 3.0])
!!v3 = vector([1.0d0, 2.0d0, 3.0d0])
!!v4 = vector(4) ! [0, 0, 0, 0]
!!v5 = vector(5, 3) ! [3, 3, 3, 3, 3]
!!v6 = vector(2, 2.0) ! [2, 2]
!!v7 = vector(val=4.7d0, dim=3) ! [4.7, 4.7, 4.7]
!!v8 = vector(v7)
!!```
    procedure :: vector_constructor_int
    procedure :: vector_constructor_r32
    procedure :: vector_constructor_r64
    procedure :: vector_constructor_dim
    procedure :: vector_constructor_dim_value_int
    procedure :: vector_constructor_dim_value_r32
    procedure :: vector_constructor_dim_value_r64
    procedure :: vector_constructor_vector

end interface 

interface operator(*)   
!! Extend multiplication operator to allow a scalar \(\alpha\) times a vector such that 
    procedure :: int_times_vector
    procedure :: r32_times_vector
    procedure :: r64_times_vector
end interface

interface operator(/)
!! Extend division operator to allow a scalar divided by a vector
    procedure :: int_div_vector
    procedure :: r32_div_vector
    procedure :: r64_div_vector
end interface

contains 


!=============================================================================!
!=                               Constructors                                =!
!=============================================================================!

    pure subroutine new_constructor(self, dim) 
    !! allocate the proper space for the elements of vector \(v\) and set the dimension to \(dim\)
        class(vector), intent(inout) :: self !! \(v\)
        integer, intent(in) :: dim !! \(n\)

        self%dim = dim
        allocate(self%v(dim))        

    end subroutine

    pure function vector_constructor_int(array) result(this)
    !! Construct a vector \(v\) from an array of integers
        integer, dimension(:), intent(in) :: array !! input data
        type(vector) :: this !! \(v\)

        call this%new_(size(array))
        this%v = array

    end function

    pure function vector_constructor_r32(array) result(this)
    !! Construct a vector \(v\) from an array of single precision reals    
        real(real32), dimension(:), intent(in) :: array !! input data
        type(vector) :: this !! \(v\)

        call this%new_(size(array))
        this%v = array

    end function

    pure function vector_constructor_r64(array) result(this)
    !! Construct a vector \(v\) from an array of double precision reals    
        real(real64), dimension(:), intent(in) :: array !! input data
        type(vector) :: this !! \(v\)

        call this%new_(size(array))
        this%v = array

    end function

    elemental function vector_constructor_dim(dim) result(this)
    !! Construct a vector by declaring its size
    !! Allocate an \(n\)-dimensional vector and fill its values with 0
        integer, intent(in) :: dim !! \(n\)
        type(vector) :: this !! \(v\)

        call this%new_(dim)
        this%v = 0

    end function

    elemental function vector_constructor_dim_value_int(dim, val) result(this)
    !! Construct a vector \(v\) of dimension \(n\) and fill its values with integer \(val\)
        integer, intent(in) :: dim !! \(n\)
        integer, intent(in) :: val !! \(val\)

        type(vector) :: this !! \(v\)

        call this%new_(dim)
        this%v = val

    end function

    elemental function vector_constructor_dim_value_r32(dim, val) result(this)
    !! Construct a vector \(v\) of dimension \(n\) and fill its values with single precision real \(val\)
        integer, intent(in) :: dim !! \(n\)
        real(real32), intent(in) :: val !! \(val\)

        type(vector) :: this !! \(v\)

        call this%new_(dim)
        this%v = val

    end function

    elemental function vector_constructor_dim_value_r64(dim, val) result(this)
    !! Construct a vector \(v\) of dimension \(n\) and fill its values with double precision real \(val\)
        integer, intent(in) :: dim !! \(n\)
        real(real64), intent(in) :: val !! \(val\)

        type(vector) :: this !! \(v\)

        call this%new_(dim)
        this%v = val

    end function

    elemental function vector_constructor_vector(v1) result(v2)
    !! Construct a vector from another vector
    !! @Note
    !!Not very efficient due do the multiple copies that occur (About three times slower than assigment)
        class(vector), intent(in) :: v1
        type(vector) :: v2

        v2 = v1

    end function

    elemental function vector_constructor_eye(self, dim, col) result(v2)
    !! Construct a vector \(v\) that is equal to the \(col\)th column of the Identity matrix \(I_{dim}\)
    !! The \(dim\) and \(col\) parameters are both optional. If the \(dim\) parameter is absent, then take the \(col\)th
    !! column from the identity matrix whose dimension is \(v_{dim}\). If the \(col\) parameter is missing, set its default
    !! value to the dimension.
    !!```fortran
    !!type(vector) :: v1, i1, i2, i3
    !!
    !!v1 = vector(4, 0) ! Initialize a vector with dimension 4 and all 0 elements
    !!
    !!i1 = v1%eye() ! Return [0, 0, 0, 1]
    !!i2 = v1%eye(col = 3) ! Return [0, 0, 1, 0]
    !!i3 = v1%eye(dim=5, col=2) ! Return [0, 1, 0, 0, 0]
    !!```
        class(vector), intent(in) :: self !! \(v\)
        integer, intent(in), optional :: dim !! Dimension of the Identity matrix
        integer, intent(in), optional :: col !! Column to extract

        integer :: dim_
        integer :: col_

        type(vector) :: v2 !! \(col\)th column of the Identity matrix \(I_{dim}\)

        if(present(dim)) then
            dim_ = dim
        else 
            dim_ = self%dim
        end if

        if(present(col)) then
            col_ = col
        else
            col_ = dim_
        end if


        v2 = vector(dim_, 0)
        call v2%set(col_, 1)          


    end function

    elemental subroutine clear_vector(self)
    !! Deallocate a vector \(v\) if it is allocated, set the dimension equal to 0
        class(vector), intent(inout) :: self !! \(v\)

        if (self%allocated()) then
            call self%dealloc_()
        else
            self%dim = 0
        end if

    end subroutine

    elemental subroutine allocate_vector_data(self, dim) 
    !! Allocate the underlying array containing \(v\)'s data and set \(v\)'s dimension to \(dim\)
        class(vector), intent(inout) :: self !! \(v\)
        integer, intent(in) :: dim !! \(n\)

        integer :: ierr

        allocate(self%v(dim), STAT=ierr)
        self%dim = dim

        if (ierr /= 0) error stop "Error allocating vector"

    end subroutine

    elemental subroutine deallocate_vector_data(self) 
    !! Deallocate the underlying array containing \(v\)'s elements
        class(vector), intent(inout) :: self !! \(v\)
        integer :: ierr

        deallocate(self%v, STAT=ierr)
        self%dim = 0

        if (ierr /= 0) error stop "Error allocating vector"

    end subroutine

!=============================================================================!
!=                         Assigment Functions                               =!
!=============================================================================!
    pure subroutine vector_from_int(self, val) 
    !! Assign a vector \(v\) to an int value. If \(v\) is already allocated, fill the elements with \(val\). If \(v\)
    !! is not already allocated, create a new vector of dimension 1 and set the element equal to \(val\)
        class(vector), intent(inout) :: self !! \(v\)
        integer, intent(in) :: val !! Value used to fill \(v\)

        if(self%allocated()) then
            self%v = val ! Copy the contents of array into self
        else
            self = vector(1, val)
        end if

    end subroutine

    pure subroutine vector_from_r32(self, val) 
    !! Assign a vector \(v\) to single precision value. If \(v\) is already allocated, fill the elements with \(val\). If \(v\)
    !! is not already allocated, create a new vector of dimension 1 and set the element equal to \(val\)
        class(vector), intent(inout) :: self !! \(v\)
        real(real32), intent(in) :: val !! Value used to fill \(v\)

        if(self%allocated()) then
            self%v = val ! Copy the contents of array into self
        else
            self = vector(1, val)
        end if

    end subroutine

    pure subroutine vector_from_r64(self, val) 
    !! Assign a vector \(v\) to a double precision value. If \(v\) is already allocated, fill the elements with \(val\). If \(v\)
    !! is not already allocated, create a new vector of dimension 1 and set the element equal to \(val\)
        class(vector), intent(inout) :: self !! \(v\)
        real(real64), intent(in) :: val !! Value used to fill \(v\)

        if(self%allocated()) then
            self%v = val ! Copy the contents of array into self
        else
            self = vector(1, val)
        end if

    end subroutine

    pure subroutine vector_from_array_int(self, array) 
    !! Assign a vector \(v\) to an array of integers. 
    !! @Note
    !! The underlying data is stored with double precision floating values, but a vector can be created from 
    !! any numeric type
        class(vector), intent(inout) :: self
        integer, dimension(:), intent(in) :: array

        self%dim = size(array)
        self%v = array ! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_array_r32(self, array) 
    !! Assign a vector \(v\) to an array of integers. 
        class(vector), intent(inout) :: self
        real(real32), dimension(:), intent(in) :: array

        self%dim = size(array)
        self%v = array! Copy the contents of array into self

    end subroutine

    pure subroutine vector_from_array_r64(self, array) 
    !! Assign a vector \(v\) to an array of integers. 
        class(vector), intent(inout) :: self
        real(real64), dimension(:), intent(in) :: array

        self%dim = size(array)
        self%v = array ! Copy the contents of array into self

    end subroutine

    elemental subroutine vector_from_vector(self, v1)
    !! Copy the elements of a vector \(v_1\) into \(v\)
        class(vector), intent(inout) :: self !! \(v\)
        class(vector), intent(in) :: v1 !! \(v_1\)

        self%dim = v1%dim
        self%v = v1%v ! Copy the array contents

    end subroutine

!=============================================================================!
!=                            State Functions                                =!
!=============================================================================!

    elemental function vector_conform(self, v2) result(bool)
    !! Check if two vectors are conforming (have the same dimension)
        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        logical :: bool !! True when self%dim == v2%dim

        bool = (self%dim == v2%dim)

    end function

    elemental function vector_is_orthogonal(self, v2, eps) result(bool)
    !! Check if two vectors are orthogonal within a certain tolerance.
        class(vector), intent(in) :: self !! \(v\)
        class(vector), intent(in) :: v2 !! \(v_2\)
        real(real64), intent(in), optional :: eps !! \(\epsilon\)

        logical :: bool !! True when \(\abs{\langle v, v2 \rangle} < \epsilon \)

        if(present(eps)) then

            bool = abs(self .dot. v2) < eps

        else

            bool = abs(self .dot. v2) < vector_epsilon

        end if

    end function

    elemental function vector_is_normal(self) result(bool)
    !! Check if a vector is normal. A vector is normal if it's length is equal to 1 (within a tolerance)
        class(vector), intent(in) :: self
        logical :: bool

        if (self%length() - 1 < vector_epsilon) then
            bool = .true.
        else 
            bool = .false.
        end if

    end function

    elemental function vector_is_allocated(self) result(bool)
    !! Check if a vector is allocated
        class(vector), intent(in) :: self

        logical :: bool

        if (allocated(self%v)) then
            bool = .true.
        else
            bool = .false.
        end if

    end function

    elemental function vector_size(self) result(n)
    !! Return the number of elements allocated for \(v\)
        class(vector), intent(in) :: self !! \(v\)
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
    
    elemental subroutine vector_normalize(self)
    !! Normalize a vector such that its euclidian norm is 1

        class(vector), intent(inout) :: self     

        real(real64) :: norm

        norm = self%length()

        self%v = self%v / norm

    end subroutine

    elemental subroutine vector_orthogonalize(self, v2)
    !! Orthogonalize a vector with respect to another
    
        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2

        type(vector) :: self_copy

        self_copy = self

        call self%proj(v2)
        call self%minus(self_copy)    

    end subroutine

    elemental subroutine vector_orthonormalize(self, v2)
    !! Orthogonalize a vector with respect to another
    
        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2

        call self%orthogonalize(v2)
        call self%normalize()  

    end subroutine

    elemental function vector_normalized(self) result(normalized_vector)
    !! Normalize a vector such that its euclidian norm is 1

        class(vector), intent(in) :: self        
        type(vector) :: normalized_vector

        real(real64) :: norm

        norm = self%length()

        normalized_vector = (self/norm)

    end function

    elemental function vector_orthogonalized(self, v2) result(v3)
    !! Orthogonalize a vector with respect to another
    
        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        v3 = self .proj. v2
        v3 = self - v3

    end function

    elemental function vector_orthonormalized(self, v2) result(v3)
    !! Orthogonalize a vector with respect to another
    
        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        v3 = self%orthogonalized(v2)
        call v3%normalize()

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

    pure function vector_outer_vector(self, v2) result(array)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2        

        real(real64), dimension(self%dim, v2%dim) :: array

    end function

    elemental function vector_householder(self, normal) result(rotated)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: normal !! MUST BE A UNIT VECTOR

        type(vector) :: rotated

        rotated = self - (2 * (self .inner. normal) * normal)

    end function

    function vector_find_householder_normal(self, destination) result(normal)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: destination

        type(vector) :: normal

        print *, "Destination set as: ", destination%data()

        normal = self - destination
        call normal%normalize()

    end function

    elemental function vector_plus_vector(self, v2) result(v3)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        if (self%conform_(v2)) then
            call v3%new_(self%dim)
            v3%v = self%v + v2%v
        else 
            error stop "Cannot add nonconforming vectors"
        end if

    end function

    elemental function vector_minus_vector(self, v2) result(v3)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        if (self%conform_(v2)) then
            call v3%new_(self%dim)
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

        call v2%new_(self%dim)

        v2%v = self%v * scalar

    end function

    elemental function vector_times_scalar_r32(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real32), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new_(self%dim)

        v2%v = self%v * scalar
    
    end function

    elemental function vector_times_scalar_r64(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real64), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new_(self%dim)

        v2%v = self%v * scalar
    
    end function

    elemental function vector_hadamard_vector(self, v2) result(v3)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        if (.not. self%conform_(v2)) error stop "Cannot multiply non_conforming vectors"

        v3 = self%v * v2%v

    end function

    elemental function vector_div_vector(self, v2) result(v3)

        class(vector), intent(in) :: self
        class(vector), intent(in) :: v2

        type(vector) :: v3

        if (.not. self%conform_(v2)) error stop "Cannot multiply non_conforming vectors"

        v3 = self%v / v2%v

    end function

    elemental function vector_div_scalar_int(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        integer, intent(in) :: scalar 
        type(vector) :: v2

        call v2%new_(self%dim)

        v2%v = self%v / scalar

    end function

    elemental function vector_div_scalar_r32(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real32), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new_(self%dim)

        v2%v = self%v / scalar
    
    end function

    elemental function vector_div_scalar_r64(self, scalar) result(v2)
    !! Multiply a vector times an integer scalar

        class(vector), intent(in) :: self
        real(real64), intent(in) :: scalar 
        type(vector) :: v2

        call v2%new_(self%dim)

        v2%v = self%v / scalar
    
    end function

    elemental function int_times_vector(scalar, vec) result(v2)

        integer, intent(in) :: scalar
        class(vector), intent(in) :: vec

        type(vector) :: v2

        v2 = vec * scalar

    end function

    elemental function r32_times_vector(scalar, vec) result(v2)

        real(real32), intent(in) :: scalar
        class(vector), intent(in) :: vec

        type(vector) :: v2

        v2 = vec * scalar

    end function

    elemental function r64_times_vector(scalar, vec) result(v2)

        real(real64), intent(in) :: scalar
        class(vector), intent(in) :: vec

        type(vector) :: v2

        v2 = vec * scalar

    end function

    elemental function int_div_vector(scalar, vec) result(v2)

        integer, intent(in) :: scalar
        class(vector), intent(in) :: vec

        type(vector) :: v2

        v2 = vec * (1._real64/scalar) 

    end function

    elemental function r32_div_vector(scalar, vec) result(v2)

        real(real32), intent(in) :: scalar
        class(vector), intent(in) :: vec

        type(vector) :: v2

        v2 = vec * (1._real64/scalar) 

    end function

    elemental function r64_div_vector(scalar, vec) result(v2)

        real(real64), intent(in) :: scalar
        class(vector), intent(in) :: vec

        type(vector) :: v2

        v2 = vec * (1._real64/scalar) 

    end function

    elemental function vector_unary_minus(self) result(v2)

        class(vector), intent(in) :: self

        type(vector) ::  v2

        v2 = self * (-1)

    end function

!=============================================================================!
!=                         Operator Subroutines                              =!
!=============================================================================!



    elemental subroutine vector_proj_vector_sub(self, v2) 
    !! Project vector self ONTO v2 \(proj_v2(self)\)
        
        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2

        real(real64) :: scalar

        scalar = (self .dot. v2) / (v2 .dot. v2)       

        self = v2 * scalar
        ! Allocate the space for a new vector
        
    end subroutine

    elemental subroutine vector_plus_vector_sub(self, v2) 

        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2


        if (self%conform_(v2)) then
            self%v = self%v + v2%v
        else 
            error stop "Cannot add nonconforming vectors"
        end if

    end subroutine

    elemental subroutine vector_householder_sub(self, normal) 

        class(vector), intent(inout) :: self
        class(vector), intent(in) :: normal !! MUST BE A UNIT VECTOR

        call self%minus(2 * (self .proj. normal))

    end subroutine

    elemental subroutine vector_minus_vector_sub(self, v2) 

        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2


        if (self%conform_(v2)) then
            self%v = self%v - v2%v
        else 
            error stop "Cannot add nonconforming vectors"
        end if

    end subroutine

    elemental subroutine vector_times_scalar_int_sub(self, scalar) 
    !! Multiply a vector times an integer scalar

        class(vector), intent(inout) :: self
        integer, intent(in) :: scalar 

        self%v = self%v * scalar

    end subroutine

    elemental subroutine vector_times_scalar_r32_sub(self, scalar) 
    !! Multiply a vector times an integer scalar

        class(vector), intent(inout) :: self
        real(real32), intent(in) :: scalar 

        self%v = self%v * scalar
    
    end subroutine

    elemental subroutine vector_times_scalar_r64_sub(self, scalar) 
    !! Multiply a vector times an integer scalar

        class(vector), intent(inout) :: self
        real(real64), intent(in) :: scalar 

        self%v = self%v * scalar
    
    end subroutine

    elemental subroutine vector_div_scalar_int_sub(self, scalar) 
    !! Multiply a vector times an integer scalar

        class(vector), intent(inout) :: self
        integer, intent(in) :: scalar 

        self%v = self%v / scalar

    end subroutine

    elemental subroutine vector_div_scalar_r32_sub(self, scalar) 
    !! Multiply a vector times an integer scalar

        class(vector), intent(inout) :: self
        real(real32), intent(in) :: scalar 

        self%v = self%v / scalar
    
    end subroutine

    elemental subroutine vector_div_scalar_r64_sub(self, scalar) 
    !! Multiply a vector times an integer scalar

        class(vector), intent(inout) :: self
        real(real64), intent(in) :: scalar 

        self%v = self%v / scalar
    
    end subroutine

    elemental subroutine vector_unary_minus_sub(self) 

        class(vector), intent(inout) :: self

       self%v = -self%v

    end subroutine

    elemental subroutine vector_times_vector_sub(self, v2)

        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2

        if (.not. self%conform_(v2)) error stop "Cannot multiply non_conforming vectors"

        self%v = self%v * v2%v

    end subroutine

    elemental subroutine vector_div_vector_sub(self, v2)

        class(vector), intent(inout) :: self
        class(vector), intent(in) :: v2

        if (.not. self%conform_(v2)) error stop "Cannot multiply non_conforming vectors"

        self%v = self%v / v2%v

    end subroutine

!=============================================================================!
!=                           Static Functions                                =!
!=============================================================================!

    elemental subroutine vector_zero(self, dim)

        class(vector), intent(inout) :: self
        integer, intent(in), optional :: dim

        if (present(dim)) then
            call self%new_(dim)
        else
            call self%new_(1)
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

        ortho(1) = vector_array(1)%normalized()

        do i = 2,k

            ortho(i) = vector_array(i)

            do j = 1,i-1

                ortho(i) = ortho(i) - (ortho(i) .proj. ortho(j))

            end do

            call ortho(i)%normalize()

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