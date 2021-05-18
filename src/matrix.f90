module matrix_m

!! A matrix can be used to mathematically represent many different structures. 
!! A matrix is a rank 1 array of vectors

use vector_m
use iso_fortran_env, only: real64, real32

implicit none
private


type, public :: matrix

    private

    integer :: k = 1 !! The number of vectors
    integer :: n = 1 !! Dimension of vectors
    type(vector), dimension(:), pointer :: m !! The vectors stored in a matrix
    logical :: m_allocated = .false.

contains 

    private

    procedure, public :: new => new_matrix
    procedure, public :: clear => clear_matrix
    procedure, public :: print => print_matrix
    procedure, public :: vec => access_vector_matrix !! Get the kth vector in the matrix
    procedure, public :: at => at_index_matrix
    procedure, public :: gram_schmidt => gram_schmidt_matrix
    procedure, public :: is_orthonormal => is_orthonormal_matrix
    
    generic, public :: set => set_int_, set_r32_, set_r64_
    generic, public :: assignment(=) => from_array_int_, from_array_r32_, from_array_r64_, from_matrix

    procedure :: from_array_int_ => matrix_from_rank2_array_int
    procedure :: from_array_r32_ => matrix_from_rank2_array_r32
    procedure :: from_array_r64_ => matrix_from_rank2_array_r64
    procedure :: from_matrix => matrix_from_matrix

    procedure :: set_int_ => set_index_matrix_int
    procedure :: set_r32_ => set_index_matrix_r32
    procedure :: set_r64_ => set_index_matrix_r64
    
    procedure :: alloc_ => allocate_matrix_data
    procedure :: dealloc_ => deallocate_matrix_data

    


end type


contains

    subroutine new_matrix(self, n, k)

        class(matrix), intent(inout) :: self
        integer, intent(in) :: n !! The dimension of each constituent vector
        integer, intent(in) :: k !! The number of vectors

        call self%clear()
 
        if(k <= 0 .or. n <= 0) error stop "Cannot instantiate vector with 0 or negative dimension"

        self%k = k
        self%n = n 

        call self%alloc_()  

    end subroutine

    subroutine matrix_from_rank2_array_int(self, array) 

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

    subroutine matrix_from_rank2_array_r32(self, array) 

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

    subroutine matrix_from_rank2_array_r64(self, array) 

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

    subroutine matrix_from_matrix(self, m) 

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

    subroutine clear_matrix(self)

        class(matrix), intent(inout) :: self

        self%k = 0
        self%n = 0

        if (self%m_allocated) then

            call self%dealloc_()

        end if

    end subroutine

    subroutine allocate_matrix_data(self) 

        class(matrix), intent(inout) :: self
        integer :: ierr

        allocate(self%m(self%k), STAT=ierr)

        call self%m%zero(self%n)

        if (ierr /= 0) error stop "Error allocating vector"

        self%m_allocated = .true.

    end subroutine

    subroutine deallocate_matrix_data(self) 

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

    function at_index_matrix(self, i, j) result(element)

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
        integer :: x
        
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

    function access_vector_matrix(self, v) result(vec)
    !! Get a copy of the vth vector

        class(matrix), intent(in) :: self
        integer, intent(in) :: v

        type(vector) :: vec

        if (v < 1 .or. v > self%k) error stop "Out of bounds index"

        vec = self%m(v)

    end function

    function gram_schmidt_matrix(self) result(ortho)

        class(matrix), intent(in) :: self
        type(matrix) :: ortho

        integer i, j

        call ortho%new(self%n, self%k)

        ! print *, "dimension of input basis = ", self%n
        ! print *, "number of basis vectors = ", self%k
        ! print *, "orthonormal_basis set to 0"

        ortho%m(1) = self%m(1)%normalize()

        do i = 2,self%k

            ortho%m(i) = self%m(i)

            do j = 1,i-1

                ortho%m(i) = ortho%m(i) - (ortho%m(i) .proj. ortho%m(j))
                ! ortho%m(i) = (ortho%m(i) .proj. ortho%m(j)) - ortho%m(i)

            end do

            ortho%m(i) = ortho%m(i)%normalize()

        end do
        
    end function
    
    function is_orthonormal_matrix(self) result(bool)

        class(matrix), intent(in) :: self
        logical :: bool

        integer :: k, i, j

    
        if(all(self%m%is_normal())) then

            do i = 1,k
                do j = i+1,k
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


end module