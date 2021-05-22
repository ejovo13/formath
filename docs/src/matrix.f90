module matrix_m
!! A matrix is a wrapper class for a rank 1 array of vector objects. <br>
!! Matrices can be used to represent a variety of mathematical structures. This class is primarily used
!! to bind a selection of Linear Algebra algorithms to a matrix object. Matrices can be instantiated by assignment 
!! of a rank2 array of any type, but the underlying data will be stored with double precision. <br>
!!
!!```fortran
!!type(matrix) :: m, ortho_basis
!!
!!m = reshape([1, 2, 3, 4], [2, 2]) ! Create a 2x2 matrix
!!print*, "M: "
!!call m%print()
!!
!!ortho_basis = m%gram_schmidt() ! Compute an orthonormal basis using the Gram-Schmidt method
!!
!!print"(A)", "Ortho:"
!!call ortho_basis%print()!!
!!```
!!
!!output:
!!
!!![test](../media/matrix_m_output.png)
!!

use vector_m
use iso_fortran_env, only: real64, real32

implicit none
private

type, public :: matrix

    private

    integer :: k = 1 !! Number of vectors
    integer :: n = 1 !! Dimension of vectors
    type(vector), dimension(:), pointer :: m !! The vectors stored in a matrix
    logical :: m_allocated = .false. !! Allocation status of pointer

contains 

    private

    procedure, public :: new => new_matrix !! Create a new matrix
    procedure, public :: clear => clear_matrix !! Clear all of the elements of a matrix 
    procedure, public :: print => print_matrix !! Print the contents of a matrix
    procedure, public :: vec => access_vector_matrix !! Get the kth vector in the matrix
    procedure, public :: at => at_index_matrix !! Get the element at the index (i, j)
    procedure, public :: gram_schmidt => gram_schmidt_matrix !! Compute an otrthonormal basis for the vector space spanned by the columns of a matrix
    procedure, public :: is_orthonormal => is_orthonormal_matrix !! Check whether a matrix is orthonormal
    procedure, public :: as_array => matrix_as_array !! Return a rank2 Fortran array
    
    generic, public :: set => set_int_, set_r32_, set_r64_ !! Set the value of \(a_{i,j})
    generic, public :: assignment(=) => from_array_int_, from_array_r32_, from_array_r64_, from_matrix !! Assign the contents of a matrix from a rank2 Fortran array

    procedure :: from_array_int_ => matrix_from_rank2_array_int
    procedure :: from_array_r32_ => matrix_from_rank2_array_r32
    procedure :: from_array_r64_ => matrix_from_rank2_array_r64
    procedure :: from_matrix => matrix_from_matrix

    procedure :: set_int_ => set_index_matrix_int
    procedure :: set_r32_ => set_index_matrix_r32
    procedure :: set_r64_ => set_index_matrix_r64
    
    procedure :: alloc_ => allocate_matrix_data !! Allocate the space for an array containing the matrix's elements
    procedure :: dealloc_ => deallocate_matrix_data  !! Deallocate the underlying container for a matrix's elements

end type


contains

    elemental subroutine new_matrix(self, n, k)
    !! Wipe the contents of a matrix and allocate the proper amount of space
        class(matrix), intent(inout) :: self !! Matrix object to wipe
        integer, intent(in) :: n !! Dimension of each constituent vector
        integer, intent(in) :: k !! Number of vectors

        call self%clear()
 
        if(k <= 0 .or. n <= 0) error stop "Cannot instantiate vector with 0 or negative dimension"

        self%k = k
        self%n = n 

        call self%alloc_()  

    end subroutine

    pure subroutine matrix_from_rank2_array_int(self, array) 
    !! Assign a matrix from a rank2 integer array
        class(matrix), intent(inout) :: self
        integer, dimension(:,:), intent(in) :: array

        integer :: i, k, n

        n = size(array, dim=1)
        k = size(array, dim=2)

        call self%new(n, k)

        do i = 1, k

            self%m(i) = array(:,i)

        end do

    end subroutine

    pure subroutine matrix_from_rank2_array_r32(self, array) 
    !! Assign a matrix from a rank2 single precision real array
        class(matrix), intent(inout) :: self
        real(real32), dimension(:,:), intent(in) :: array

        integer :: i, k, n

        n = size(array, dim=1)
        k = size(array, dim=2)

        call self%new(n, k)

        do i = 1, k

            self%m(i) = array(:,i)

        end do

    end subroutine

    pure subroutine matrix_from_rank2_array_r64(self, array) 
    !! Assign a matrix from a rank2 double precision array
        class(matrix), intent(inout) :: self
        real(real64), dimension(:,:), intent(in) :: array

        integer :: i, k, n

        n = size(array, dim=1)
        k = size(array, dim=2)

        call self%new(n, k)

        do i = 1, k

            self%m(i) = array(:,i)

        end do

    end subroutine

    elemental subroutine matrix_from_matrix(self, m) 
    !! 
        class(matrix), intent(inout) :: self
        class(matrix), intent(in) :: m

        integer :: i

        self%n = m%n
        self%k = m%k

        call self%new(m%n, m%k)

        do i = 1, self%k

            self%m(i) = m%vec(i)

        end do

    end subroutine

    elemental subroutine clear_matrix(self)

        class(matrix), intent(inout) :: self

        self%k = 0
        self%n = 0

        if (self%m_allocated) then

            call self%dealloc_()

        end if

    end subroutine

    elemental subroutine allocate_matrix_data(self) 

        class(matrix), intent(inout) :: self
        integer :: ierr

        allocate(self%m(self%k), STAT=ierr)

        call self%m%zero(self%n)

        if (ierr /= 0) error stop "Error allocating vector"

        self%m_allocated = .true.

    end subroutine

    elemental subroutine deallocate_matrix_data(self) 

        class(matrix), intent(inout) :: self
        integer :: ierr

        call self%m%dealloc_()

        deallocate(self%m, STAT=ierr)

        if (ierr /= 0) error stop "Error allocating vector"

        self%m_allocated = .false.

    end subroutine

    subroutine print_matrix(self)

        class(matrix), intent(in) :: self

        character(10) :: num_fmt
        character(100) :: fmt

        integer :: i, j

        write(num_fmt,"(I0)") self%k


        fmt = "(" // num_fmt // "(G11.5, 2X))"
        ! fmt = "(" // num_fmt // "(G0, 2X))"

        do i = 1, self%n

            print fmt, (self%at(i, j), j = 1, self%k)

        end do

    end subroutine

    elemental function at_index_matrix(self, i, j) result(element)

        class(matrix), intent(in) :: self
        integer, intent(in) :: i !! ith element
        integer, intent(in) :: j !! jth vector
        
        real(real64) :: element
        
        element = self%m(j)%at(i)

    end function

    subroutine set_index_matrix_int(self, i, j, x) 

        class(matrix), intent(in) :: self
        integer, intent(in) :: i !! ith element
        integer, intent(in) :: j !! jth vector        
        integer, intent(in) :: x
        
        call self%m(j)%set(i, x)

    end subroutine

    subroutine set_index_matrix_r32(self, i, j, x) 

        class(matrix), intent(in) :: self
        integer, intent(in) :: i !! ith element
        integer, intent(in) :: j !! jth vector        
        real(real32) :: x
        
        call self%m(j)%set(i, x)

    end subroutine

    subroutine set_index_matrix_r64(self, i, j, x) 

        class(matrix), intent(in) :: self
        integer, intent(in) :: i !! ith element
        integer, intent(in) :: j !! jth vector        
        real(real64) :: x
        
        call self%m(j)%set(i, x)

    end subroutine

    elemental function access_vector_matrix(self, v) result(vec)
    !! Get a copy of the vth vector

        class(matrix), intent(in) :: self
        integer, intent(in) :: v

        type(vector) :: vec

        if (v < 1 .or. v > self%k) error stop "Out of bounds index"

        vec = self%m(v)

    end function

    elemental function gram_schmidt_matrix(self) result(ortho)

        class(matrix), intent(in) :: self
        type(matrix) :: ortho

        integer i, j, n, k

        n = self%n
        k = self%k 

        if(k > n) k = n !! If there are more vectors than the dimension of the vector, only output n vectors

        call ortho%new(n, k)

        ! print *, "dimension of input basis = ", self%n
        ! print *, "number of basis vectors = ", self%k
        ! print *, "orthonormal_basis set to 0"

        ortho%m(1) = self%m(1)%normalized()

        do i = 2,k

            ortho%m(i) = self%m(i)

            do j = 1,i-1

                ortho%m(i) = ortho%m(i) - (ortho%m(i) .proj. ortho%m(j))
                ! ortho%m(i) = (ortho%m(i) .proj. ortho%m(j)) - ortho%m(i)

            end do

            call ortho%m(i)%normalize()

        end do
        
    end function
    
    elemental function is_orthonormal_matrix(self) result(bool)

        class(matrix), intent(in) :: self
        logical :: bool

        integer :: i, j

    
        if(all(self%m%is_normal())) then

            do i = 1,self%k
                do j = i+1,self%k
                    if (.not. self%m(i)%is_ortho(self%m(j))) then
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

    pure function matrix_as_array(self) result(array)

        class(matrix), intent(in) :: self
        real(real64), dimension(self%n, self%k) :: array

        integer :: i

        forall(i = 1:self%k) 

            array(i,:) = self%m(i)%data()

        end forall

    end function


end module