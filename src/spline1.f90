subroutine spline1(x,f,fs,n)
        use mod_constantes, only : dp
  implicit real(kind=dp) (a-h,o-z)
!integer, parameter :: dp = selected_real_kind(15, 307)
real(kind=dp) , allocatable, dimension(:,:) :: a
real(kind=dp) , intent(in), dimension(0:n)   :: x, f
real(kind=dp) , intent(out), dimension(0:n)   :: fs
real(kind=dp) , allocatable, dimension(:)   :: b, h
integer  ,  intent(in)     :: n
integer        :: i,j,k
real(kind=dp) :: z,xp
allocate(a(0:n,1:3))
allocate(b(0:n))
allocate(h(0:n))
a(:,:)=0.d0
do i=0,n-1
  h(i)=x(i+1)-x(i)
enddo
do i=1,n-1
  a(i,2)=2.d0*(h(i)+h(i-1))
enddo
a(0,2)=2.d0*h(0)
a(n,2)=2.d0*h(n-1)
do i=0,n-1
  a(i+1,1)=h(i)
enddo
!a(n,1)=0.d0
do i=1,n
  a(i-1,3)=h(i-1)
enddo
!a(0,3)=0.d0
do i=1,n-1
  b(i)=6.d0*((f(i+1)-f(i))/h(i)-(f(i)-f(i-1))/h(i-1))
enddo
  b(0)=6.d0*((f(1)-f(0))/h(0)-((f(i)-f(i-1))/h(i-1)+((f(1)-f(0))/h(0)-(f(2)-f(1))/h(1))))
  b(n)=-6.d0*((f(n)-f(n-1))/h(n-1))

   if(f(n-1)*f(n)<=0.d0.or.f(n-1)/f(n)<=1.0d0.or.dabs(f(n))<1.d-2)b(n)=0.d0
!   write(9,*)b(0),b(n)
!b(0)=0.d0
!b(n)=0.d0
do i=1,n
  a(i,1)=a(i,1)/a(i-1,2)
  a(i,2)=a(i,2)-a(i,1)*a(i-1,3)
enddo 
do i=1,n
  b(i)=b(i)-a(i,1)*b(i-1)
enddo
fs(n)=b(n)/a(n,2)
do i=n-1,0,-1
  fs(i)=(b(i)-a(i,3)*fs(i+1))/a(i,2)
enddo
end subroutine spline1
