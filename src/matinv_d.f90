recursive subroutine matinv_d(s,n)
 use OMP_LIB
  use cublas
  use cusolverDn
  use cudafor
integer :: i, j, k, l
integer,intent(in) :: n
real(kind=8), intent(inout), dimension(n,n) :: s

real(kind=8), device, allocatable, dimension(:,:) :: a_d, b_d, d_d, e_d, f_d, h_d
real(kind=8), allocatable, dimension(:,:) :: a,b
  integer,allocatable,dimension(:)::isuppz
  real(kind=8),allocatable,dimension(:)::WORK
  logical iprint
iprint=.true.
imp=9
!n=size(a,1)
  nn=n/2
  m=n-nn
  LWORK=64*m
  allocate(WORK(LWORK))
  allocate(isuppz(m))
  allocate(a_d(nn,nn),e_d(nn,nn),b_d(nn,m),f_d(nn,m),d_d(m,m),h_d(m,m))
  allocate(a(nn,nn),b(m,m))
!sm1(:,:)=0.0d00
!
! 
!
!
!
!
!
   b_d=s(1:nn,nn+1:n)
!
! b_d  ------>  B
!
   d_d=s(nn+1:n,nn+1:n)
!
! d_d  ------>  D
!
!
!
!
!
!
   a_d=s(1:nn,1:nn)
!   write(imp,*)s_d(i,j)
!write(*,*)'s'
!write(*,*) s_d
!write(*,*)
!          call imprime2(iprint,imp,nn,a_d,'a_d       ')
!          flush(imp)
!
! a_d  ------>  A
!
a=a_d
if(nn<500)then
     call dsytrf('U',nn,a,nn,isuppz,work,lwork,info)
     call dsytri('U',nn,a,nn,isuppz,work,info)
else
call matinv(a,nn)
endif
a_d=a
!do j=1,nn
!!$OMP PARALLEL DO PRIVATE (i,j)
!do i=j,nn
!   a_d(i,j)=a_d(j,i)
!enddo
!!$OMP END PARALLEL DO
!enddo
!
! a_d  ------>  A-1                                 V
!
!
!
!
     call dsymm('L','U',nn,m,1.d00,a_d,nn,b_d,nn,0.d0,f_d,nn)
!     call ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
!
! f_d  ------>  A-1 * B                             X
!
     call dgemm('T','N',m,m,nn,-1.d00,b_d,nn,f_d,nn,0.d0,h_d,m)
!$cuf kernel do(2) <<<*,*>>>
do i=1,m
do j=1,m
   h_d(i,j)=d_d(i,j)+h_d(i,j)
enddo
enddo
!
! h_d  ------>  D - C * A-1 * B                     
!
!          call imprime2(iprint,imp,m,h_d,'h_d       ')
!h_d=a_d
b=h_d
if(m<500)then
     call dsytrf('U',m,b,m,isuppz,work,lwork,info)
     call dsytri('U',m,b,m,isuppz,work,info)
else
call matinv(b,m)
endif
h_d=b
!          call imprime2(iprint,imp,m,h_d,'h_d       ')
!!$cuf kernel do(2) <<<*,*>>>
!do j=1,m
!do i=j,m
!   h_d(i,j)=h_d(j,i)
!enddo
!enddo
!
! 
!
!
!
   s(nn+1:n,nn+1:n)=h_d
!
! h_d  ------>  (D - C * A-1 * B)-1  =  H
!
!
!
!
!
     call dsymm('R','U',nn,m,-1.d00,h_d,m,f_d,nn,0.d0,b_d,nn)
!
!
!
!
!
!
   s(1:nn,nn+1:n)=b_d
!$OMP PARALLEL DO PRIVATE (i,j)
do j=1,m
do i=1,nn
   s(j+nn,i)=s(i,j+nn)
enddo
enddo
!$OMP END PARALLEL DO
!
! b_d  ------> - A-1 * B * (D - C * A-1 * B)-1  =  F                  -XH
!
!
!
!
     call dgemm('N','T',nn,nn,m,-1.d00,b_d,nn,f_d,nn,0.d0,e_d,nn)
!
!
!
!$cuf kernel do(2) <<<*,*>>>
do i=1,nn
do j=1,nn
   e_d(i,j)=a_d(i,j)+e_d(i,j)
enddo
enddo
   s(1:nn,1:nn)=e_d
!$OMP PARALLEL DO PRIVATE (i,j)
do j=1,nn
do i=1,j
   s(j,i)=s(i,j)
enddo
enddo
!$OMP END PARALLEL DO
!          call imprime2(iprint,imp,nn,e_d,'e_d       ')
!
! e_d  ------> V + Z * Xt    =   E 
!

!s=sm1
deallocate(a_d,b_d,d_d,e_d,f_d,h_d)
deallocate(isuppz,WORK,a,b)
end subroutine matinv_d
