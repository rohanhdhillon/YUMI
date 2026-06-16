module mod_utildev
contains
attributes(global) subroutine daddmat_d(a, b)
implicit none
real(8) :: a(:,:), b(:,:)
integer :: i, j
i = (blockIdx%x-1)*blockDim%x + threadIdx%x
j = (blockIdx%y-1)*blockDim%y + threadIdx%y
if (i<=size(a,1) .and. j<=size(a,2)) &
a(i,j) = a(i,j) + b(i,j)
end subroutine daddmat_d
attributes(global) subroutine dminusmat_d(a, b)
implicit none
real(8) :: a(:,:), b(:,:)
integer :: i, j
i = (blockIdx%x-1)*blockDim%x + threadIdx%x
j = (blockIdx%y-1)*blockDim%y + threadIdx%y
if (i<=size(a,1) .and. j<=size(a,2)) &
a(i,j) = -a(i,j) + b(i,j)
end subroutine dminusmat_d
attributes(global) subroutine damat_d(a, b)
implicit none
real(8) :: b(:,:)
real(8), value :: a
integer :: i, j
i = (blockIdx%x-1)*blockDim%x + threadIdx%x
j = (blockIdx%y-1)*blockDim%y + threadIdx%y
if (i<=size(b,1) .and. j<=size(b,2)) &
b(i,j) = a * b(i,j)
end subroutine damat_d
attributes(global) subroutine dunmat_d(a)
implicit none
real(8) :: a(:,:)
integer :: i, j
i = (blockIdx%x-1)*blockDim%x + threadIdx%x
j = (blockIdx%y-1)*blockDim%y + threadIdx%y
if (i<=size(a,1) .and. j<=size(a,2)) &
a(i,j) = 0.d00
i = (blockIdx%x-1)*blockDim%x + threadIdx%x
if (i<=size(a,1)) &
a(i,i) = 1.d00
end subroutine dunmat_d
end module mod_utildev
