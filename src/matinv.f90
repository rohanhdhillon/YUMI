recursive subroutine matinv(s,n)

   !!This is implementation of the algorithm from https://ieeexplore.ieee.org/document/7300068
    
   use OMP_LIB
   ! this subroutine calls lapack for smaller matrices.  No need to magma, I think. 

   integer :: i, j, k, l
   integer,intent(in) :: n
   real(kind=8), intent(inout), dimension(n,n) :: s

   real(kind=8), allocatable, dimension(:,:) :: a_d, b_d, d_d, e_d, f_d, h_d
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
   !sm1(:,:)=0.0d00
   !
   ! 
   !
   !
   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
     do i=1,nn
       b_d(i,j)=s(i,j+nn)
     enddo
   enddo
   !$OMP END PARALLEL DO
   !
   ! b_d  ------>  B
   !
   
   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
     do i=1,j
      d_d(i,j)=s(i+nn,j+nn)
     enddo
   enddo
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
     do i=1,j
       d_d(j,i)=d_d(i,j)
     enddo
   enddo
   !$OMP END PARALLEL DO

   !
   ! d_d  ------>  D
   !
   !
   !
   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,nn
     do i=1,j
       a_d(i,j)=s(i,j)
     enddo
   enddo
   !$OMP END PARALLEL DO

  
   !
   ! a_d  ------>  A
   !
   if(nn<500)then
      call dsytrf('U',nn,a_d,nn,isuppz,work,lwork,info)
      call dsytri('U',nn,a_d,nn,isuppz,work,info)
      
      !Following does not work for some reason
      !call dgetrf(nn,nn,a_d,nn,isuppz,info)
      !call dgetri(nn,a_d,nn,isuppz,work,lwork,info)
   else  
      call matinv(a_d,nn)
   endif
   

   !If you use dgetrf you don't need this one
   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,nn
     do i=1,j
       a_d(j,i)=a_d(i,j)
     enddo
   enddo
   !$OMP END PARALLEL DO
   !
   ! a_d  ------>  A-1                                 V
  

   !
     call dsymm('L','U',nn,m,1.d00,a_d,nn,b_d,nn,0.d0,f_d,nn)
   !     call ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
   !
   ! f_d  ------>  A-1 * B                             X
   !
   call dgemm('T','N',m,m,nn,-1.d00,b_d,nn,f_d,nn,0.d0,h_d,m)

   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
      do i=1,m
        h_d(i,j)=d_d(i,j)+h_d(i,j)
      enddo
   enddo
   !$OMP END PARALLEL DO
  
   !h_d=a_d

   if(m<500)then
     call dsytrf('U',m,h_d,m,isuppz,work,lwork,info)
     call dsytri('U',m,h_d,m,isuppz,work,info)
     
     !Following does not work for some reason!
     !call dgetrf(m, m,h_d,m,isuppz,info)
     !call dgetri(m, h_d, m,isuppz,work,lwork,info)
     

   else
     call matinv(h_d,m)
   endif

   !          call imprime2(iprint,imp,m,h_d,'h_d       ')

   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
      do i=1,j
       h_d(j,i)=h_d(i,j)
     enddo
   enddo
   !$OMP END PARALLEL DO

   !
   ! 
   !
   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
     do i=1,m
       s(i+nn,j+nn)=h_d(i,j)
     enddo
   enddo

   !$OMP END PARALLEL DO
   !
   ! h_d  ------>  (D - C * A-1 * B)-1  =  H
   !
   !
   call dsymm('R','U',nn,m,-1.d00,h_d,m,f_d,nn,0.d0,b_d,nn)
   !
   !
   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,m
     do i=1,nn
       s(i,j+nn)=b_d(i,j)
       !   s(j+nn,i)=b_d(i,j)
     enddo
   enddo

   !$OMP END PARALLEL DO

   !
   ! b_d  ------> - A-1 * B * (D - C * A-1 * B)-1  =  F                  -XH
   !
   !
   !
   call dgemm('N','T',nn,nn,m,-1.d00,b_d,nn,f_d,nn,0.d0,e_d,nn)
   !
   !
   !$OMP PARALLEL DO PRIVATE (i,j)
   do i=1,nn
     do j=1,nn
       e_d(i,j)=a_d(i,j)+e_d(i,j)
     enddo
   enddo
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO PRIVATE (i,j)
   do j=1,nn
     do i=1,j
       s(i,j)=e_d(i,j)
       ! s(j,i)=e_d(i,j)
     enddo
   enddo  
   !$OMP END PARALLEL DO

   !          call imprime2(iprint,imp,nn,e_d,'e_d       ')
   !
   ! e_d  ------> V + Z * Xt    =   E 
   !

   !s=sm1
   deallocate(a_d,b_d,d_d,e_d,f_d,h_d)
   deallocate(isuppz,WORK)
end subroutine matinv
