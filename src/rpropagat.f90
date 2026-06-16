subroutine rpropagat(rmin,rmax,RP,nbv,n,fl)
        use mod_constantes
        use mod_pot
        use mod_base, only : ltab, jtab, k2tab, lamax, nmax
!$      use OMP_LIB        
implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,il,iu,M,l,lda,ldz, lam
  real(kind=dp),allocatable,dimension(:,:)::w,T,Tm1
  real(kind=dp),allocatable,dimension(:,:)::r1,r2,gg,gg1
  real(kind=dp),intent(out),dimension(nbv,nbv)::RP
  real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
  integer,allocatable,dimension(:)::iwork,isuppz
  real(kind=dp),allocatable,dimension(:)::WORK,lambda,TT
  real(kind=dp)::mu,E,k,rm,r0,h,h2,h1,h3,abstol,vl,vu,norm,lami
  integer , intent(in) :: nbv, n
  real(kind=dp), intent(in)::rmin,rmax
  real::result
  real,dimension(2)::tarray
  logical :: iprint
  iprint=.false.
  lda=nbv
  ldz=nbv
  !LWORK=-1
  !LIWORK=-1
  M=0
  LWORK=64*nbv
  LIWORK=64*nbv
  abstol=1.d-06
  !n=400
  rm=rmax
  r0=rmin
  h=(rm-r0)/dfloat(n)                                      ! Pas d'integration
!  call dtime(tarray, result)                                           ! Allocation de la mémoire
!        write(9,*) 'nbv =',nbv
  allocate(w(nbv,nbv))
  allocate(T(nbv,nbv))
  allocate(Tm1(nbv,nbv))
  allocate(r1(nbv,nbv))
  allocate(r2(nbv,nbv))
  allocate(gg(nbv,nbv))
  allocate(gg1(nbv,nbv))


  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  allocate(lambda(nbv))
  allocate(TT(nbv))
  allocate(isuppz(2*nbv))
!  open(12,file='matrix2.res',status='unknown')

  r=r0
r=r+0.5d0*h
!        write(9,*) r
!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
  w(i,j)=0.d0
enddo
enddo
!$OMP END PARALLEL DO
  call pot
  do lam=0,lamax
!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
     w(i,j)=w(i,j)+fl(i,j,lam)*vvl(lam)
enddo
enddo
!$OMP END PARALLEL DO
  enddo
!$OMP PARALLEL DO PRIVATE (i)
  do i=1,nbv
     w(i,i)=w(i,i)+dfloat(ltab(i)*(ltab(i)+1))/(r*r)-k2tab(i)
  enddo
!$OMP END PARALLEL DO
  call DSYEV('V','L',nbv,w,lda,lambda, work,LWORK,info) 
!  call DSYEVR('V','A','L',nbv,w,lda,vl,vu,il,iu,abstol,M,lambda, &
!       T,ldz,isuppz,work,LWORK,iwork,LiWORK,info) 
!        write(9,*)'3=', r
  T=w
  Tm1=T

!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
  RP(i,j)=0.0d00
enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE (i)
  do i=1,nbv
     RP(i,i)=1.d00/dsqrt(dabs(lambda(i)))
  enddo
!$OMP END PARALLEL DO

!  RP=matmul(RP,T)
!  T=transpose(T)
!  RP=matmul(T,RP)
call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,T,nbv,0.d00,gg,nbv)
call dgemm('T','N',nbv,nbv,nbv,1.d00,T,nbv,gg,nbv,0.d00,RP,nbv)
!RP=gg1
!r=r-0.5d0*h
!r=r-1.0d0*h
  do ir=1,n-1
     r=r+h
!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
  w(i,j)=0.d0
enddo
enddo
!$OMP END PARALLEL DO
     call pot
     do lam=0,lamax
!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
     w(i,j)=w(i,j)+fl(i,j,lam)*vvl(lam)
enddo
enddo
!$OMP END PARALLEL DO
     enddo
!$OMP PARALLEL DO PRIVATE (i)
     do i=1,nbv
        w(i,i)=w(i,i)+dfloat(ltab(i)*(ltab(i)+1))/(r*r)-k2tab(i)
     enddo
!$OMP END PARALLEL DO
!     call imprime2(iprint,9,nbv,w,'w         ')

  call DSYEV('V','L',nbv,w,lda,lambda, work,LWORK,info) 
!     call DSYEVR('V','A','L',nbv,w,lda,vl,vu,il,iu,abstol,M,lambda, &
!          T,ldz,isuppz,work,LWORK,iwork,LiWORK,info) 
        T=w
     Tm1=T*Tm1
     TT=sum(Tm1,dim=2)
     do j=1,nbv
        if(TT(j).LT.0.0d00) then
           T(:,j)=-T(:,j)
        endif
     enddo

     Tm1=T
 
     r1=0.0d00
     r2=0.0d00
     do i=1,nbv
        lami=dsqrt(dabs(lambda(i)))
        if (lambda(i).GT.0.d00) then
           r1(i,i)=1.d00/(lami*dtanh(h*lami))
           r2(i,i)=1.d00/(lami*dsinh(h*lami))
        else
           r1(i,i)=-1.d00/(lami*dtan(h*lami))
           r2(i,i)=-1.d00/(lami*dsin(h*lami))
        endif
     enddo
!     RP=matmul(RP,T)
!     T=transpose(T)
!     RP=matmul(T,RP)
!     RP=RP+r1
call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,T,nbv,0.d00,gg,nbv)
call dgemm('T','N',nbv,nbv,nbv,1.d00,T,nbv,gg,nbv,0.d00,gg1,nbv)
!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
     RP(i,j)=gg1(i,j)+r1(i,j)
enddo
enddo
!$OMP END PARALLEL DO
     call dgetrf(nbv,nbv,RP,nbv,isuppz,info)
     call dgetri(nbv,RP,nbv,isuppz,work,lwork,info)
!     RP=matmul(RP,r2)
!     RP=matmul(r2,RP)
!     RP=-RP+r1
call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,r2,nbv,0.d00,gg,nbv)
call dgemm('N','N',nbv,nbv,nbv,-1.d00,r2,nbv,gg,nbv,1.d00,r1,nbv)
!$OMP PARALLEL DO PRIVATE (i , j )
do j=1,nbv
do i=1,nbv
     RP(i,j)=r1(i,j)
enddo
enddo
!$OMP END PARALLEL DO
!     RP=matmul(T,RP)
!     T=transpose(T)
!     RP=matmul(RP,T)
call dgemm('N','N',nbv,nbv,nbv,1.d00,T,nbv,RP,nbv,0.d00,gg,nbv)
call dgemm('N','T',nbv,nbv,nbv,1.d00,gg,nbv,T,nbv,0.d00,RP,nbv)

  enddo
!call dgemm('T','N',nbv,nbv,nbv,1.d00,T,nbv,RP,nbv,0.d00,gg,nbv)
!call dgemm('N','N',nbv,nbv,nbv,1.d00,gg,nbv,T,nbv,0.d00,RP,nbv)
     call dgetrf(nbv,nbv,RP,nbv,isuppz,info)
     call dgetri(nbv,RP,nbv,isuppz,work,lwork,info)

end subroutine rpropagat
