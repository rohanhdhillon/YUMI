subroutine matinv5(s,n)
!$ use OMP_LIB
  use cublas
  use cusolverDn
  use cudafor
integer :: i, j, k, l
integer,intent(in) :: n
real(kind=8), intent(inout), dimension(n,n) :: s
!real(kind=8),  dimension(n,n) :: sm1
!real(kind=8), dimension(size(a,1),size(a,2)) :: c
!real(kind=8), dimension(n,n) :: c
!real(kind=8), intent(in), dimension(n,n) :: a,b
!real(kind=8), intent(out), dimension(n,n) :: c
!real(kind=8), dimension(size(a,1)/2,size(a,1)/2) :: s1m,s2m,s3m,s4m,t1m,t2m,t3m,t4m
!real(kind=8), dimension(size(a,1)/2,size(a,1)/2) :: p1m,p2m,p3m,p4m,p5m,p6m,p7m
!real(kind=8), dimension(size(a,1)/2,size(a,1)/2) :: u2m,u3m,u4m
!real(kind=8), dimension(size(a,1)/4,size(a,1)/4) :: s1m2,s2m2,s3m2,s4m2,t1m2,t2m2,t3m2,t4m2
!real(kind=8), dimension(size(a,1)/4,size(a,1)/4) :: p1m2,p2m2,p3m2,p4m2,p5m2,p6m2,p7m2
!real(kind=8), dimension(size(a,1)/4,size(a,1)/4) :: u2m2,u3m2,u4m2

real(kind=8), device, allocatable, dimension(:,:) :: a_d, b_d, d_d, e_d, f_d, h_d
real(kind=8), allocatable, dimension(:,:) :: a, b, d, e, f, h
!real(kind=8), dimension(n/2,n/2) :: p1m,p2m,p3m,p4m,p5m,p6m,p7m
!real(kind=8), dimension(n/2,n/2) :: u2m,u3m,u4m
!real(kind=8), dimension(n/4,n/4) :: s1m2,s2m2,s3m2,s4m2,t1m2,t2m2,t3m2,t4m2
!real(kind=8), dimension(n/4,n/4) :: p1m2,p2m2,p3m2,p4m2,p5m2,p6m2,p7m2
!real(kind=8), dimension(n/4,n/4) :: u2m2,u3m2,u4m2
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
  allocate(a(nn,nn),e(nn,nn),b(nn,m),f(nn,m),d(m,m),h(m,m))
!sm1(:,:)=0.0d00
!
! 
!
!$OMP PARALLEL DO PRIVATE (i,j)
do j=1,nn
do i=1,j
   a(i,j)=s(i,j)
!   write(imp,*)s_d(i,j)
enddo
enddo
!$OMP END PARALLEL DO
!write(*,*)'s'
!write(*,*) s_d
!write(*,*)
!          call imprime2(iprint,imp,nn,a_d,'a_d       ')
!          flush(imp)
!
! a_d  ------>  A
!
!!$OMP PARALLEL DO PRIVATE (i,j)
!do j=1,m
!do i=1,nn
!   b_d(i,j)=s(i,j+nn)
!enddo
!enddo
!!$OMP END PARALLEL DO
b_d=s(1:nn,nn+1:n)
!
! b_d  ------>  B
!
     call dsytrf('U',nn,a,nn,isuppz,work,lwork,info)
     call dsytri('U',nn,a,nn,isuppz,work,info)
do j=1,nn
!$OMP PARALLEL DO
do i=j,nn
   a(i,j)=a(j,i)
enddo
!$OMP END PARALLEL DO
enddo
!
! a_d  ------>  A-1                                 V
!
a_d=a
     call dsymm('L','U',nn,m,1.d00,a_d,nn,b_d,nn,0.d0,f_d,nn)
!     call ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
!
! f_d  ------>  A-1 * B                             X
!
!$OMP PARALLEL DO PRIVATE (i,j)
do j=1,m
do i=1,j
   d(i,j)=s(i+nn,j+nn)
enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE (i,j)
do j=1,m
do i=1,j
   d(j,i)=d(i,j)
enddo
enddo
!$OMP END PARALLEL DO
!
! d_d  ------>  D
!
d_d=d
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
h=h_d
     call dsytrf('U',m,h,m,isuppz,work,lwork,info)
     call dsytri('U',m,h,m,isuppz,work,info)
!          call imprime2(iprint,imp,m,h_d,'h_d       ')
do j=1,m
!$OMP PARALLEL DO PRIVATE (i)
do i=j,m
   h(i,j)=h(j,i)
enddo
!$OMP END PARALLEL DO
enddo
h_d=h
!$OMP PARALLEL DO PRIVATE (i,j)
do j=1,m
do i=1,m
   s(i+nn,j+nn)=h(i,j)
enddo
enddo
!$OMP END PARALLEL DO
!
! h_d  ------>  (D - C * A-1 * B)-1  =  H
!
     call dsymm('R','U',nn,m,-1.d00,h_d,m,f_d,nn,0.d0,b_d,nn)
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
     call dgemm('N','T',nn,nn,m,-1.d00,b_d,nn,f_d,nn,0.d0,e_d,nn)
!$cuf kernel do(2) <<<*,*>>>
do i=1,nn
do j=1,nn
   e_d(i,j)=a_d(i,j)+e_d(i,j)
enddo
enddo
!!$OMP PARALLEL DO PRIVATE (i,j)
!do j=1,nn
!do i=1,j
!   s(i,j)=e_d(i,j)
!   s(j,i)=e_d(i,j)
!enddo
!enddo
!!$OMP END PARALLEL DO
s(1:nn,1:nn)=e_d
!          call imprime2(iprint,imp,nn,e_d,'e_d       ')
!
! e_d  ------> V + Z * Xt    =   E 
!

!s=sm1
deallocate(a_d,b_d,d_d,e_d,f_d,h_d)
deallocate(isuppz,WORK)
end subroutine matinv5
