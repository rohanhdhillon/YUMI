!include 'lapack.f90'
!      call ijohnson(rmin,rmax,ri,fi,fs,h,fl,nmax,RP,nlev,lamax,jtab,ltab,ktab,npas)
subroutine rpropagat(rmin,rmax,ri,fi,fs,hh,fl,nmax,RP,nbv,lamax,jtab,ltab,ktab,npas)
  !  .. "Use Statements" ..
  !      USE mkl95_PRECISION, ONLY: WP => DP
  !      USE mkl95_LAPACK, ONLY: SYEV
  implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,il,iu,M,l,lda,ldz
  real(kind=8),allocatable,dimension(:,:)::w,T,Tm1
  real(kind=8),allocatable,dimension(:,:)::r1,r2,gg,gg1
  real(kind=8),intent(out),dimension(nbv,nbv)::RP
  real(kind=8),intent(in),dimension(nbv,nbv,0:lamax)::fl
  integer, intent(in), dimension(nbv) :: jtab, ltab
  real(kind=8), intent(in), dimension(nbv) :: ktab
  integer,allocatable,dimension(:)::iwork,isuppz
  real(kind=8),allocatable,dimension(:)::WORK,lambda,TT
  real(kind=8)::mu,E,k,rm,r0,r,h,h2,h1,h3,pi,abstol,vl,vu,norm,lami
  integer , intent(in) :: nbv, lamax, npas
  real(kind=8),intent(in), dimension(0:nmax) :: ri,hh
  real(kind=8),intent(in), dimension(0:nmax,0:lamax) :: fi,fs
  real(kind=8), intent(in)::rmin,rmax
  real(kind=8)::vvl
  real::result
  real,dimension(2)::tarray
  logical :: iprint
  iprint=.false.
  pi=4.d0*datan(1.d0)
  lda=nbv
  ldz=nbv
  !LWORK=-1
  !LIWORK=-1
  M=0
  LWORK=64*nbv
  LIWORK=64*nbv
  abstol=1.d-06
  n=npas
  rm=rmax
  r0=rmin
  h=(rm-r0)/dfloat(n)                                      ! Pas d'integration
  call dtime(tarray, result)                                           ! Allocation de la mémoire
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
  w(:,:)=0.d0
  do lam=0,lamax
     call pot(ri,fi(:,lam),fs(:,lam),hh,nmax,r,vvl)
     w(:,:)=w(:,:)+fl(:,:,lam)*vvl
  enddo
  do i=1,nbv
     w(i,i)=w(i,i)+dfloat(ltab(i)*(ltab(i)+1))/(r*r)-ktab(i)
  enddo
  call DSYEVR('V','A','L',nbv,w,lda,vl,vu,il,iu,abstol,M,lambda, &
       T,ldz,isuppz,work,LWORK,iwork,LiWORK,info) 
  Tm1=T

  RP(:,:)=0.0d00
  do i=1,nbv
     RP(i,i)=1.d00/dsqrt(dabs(lambda(i)))
  enddo

!  RP=matmul(RP,T)
!  T=transpose(T)
!  RP=matmul(T,RP)
call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,T,nbv,0.d00,gg,nbv)
call dgemm('T','N',nbv,nbv,nbv,1.d00,T,nbv,gg,nbv,0.d00,gg1,nbv)
RP=gg1

  do ir=1,n
     r=r+h
     w(:,:)=0.d0
     do lam=0,lamax
        call pot(ri,fi(:,lam),fs(:,lam),hh,nmax,r,vvl)
        w(:,:)=w(:,:)+fl(:,:,lam)*vvl
!        write(9,*) r,vvl
     enddo
!     call imprime2(iprint,9,nbv,w,'w         ')
!        write(16,*) r,fl(3,3,3)
     do i=1,nbv
        w(i,i)=w(i,i)+dfloat(ltab(i)*(ltab(i)+1))/(r*r)-ktab(i)
     enddo
!     call imprime2(iprint,9,nbv,w,'w         ')

     call DSYEVR('V','A','L',nbv,w,lda,vl,vu,il,iu,abstol,M,lambda, &
          T,ldz,isuppz,work,LWORK,iwork,LiWORK,info) 
     Tm1=T*Tm1
     TT=sum(Tm1,dim=1)
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
     RP=gg1+r1
     call dgetrf(nbv,nbv,RP,nbv,isuppz,info)
     call dgetri(nbv,RP,nbv,isuppz,work,lwork,info)
!     RP=matmul(RP,r2)
!     RP=matmul(r2,RP)
!     RP=-RP+r1
call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,r2,nbv,0.d00,gg,nbv)
call dgemm('N','N',nbv,nbv,nbv,-1.d00,r2,nbv,gg,nbv,0.d00,gg1,nbv)
     RP=gg1+r1
!     RP=matmul(T,RP)
!     T=transpose(T)
!     RP=matmul(RP,T)
call dgemm('N','N',nbv,nbv,nbv,1.d00,T,nbv,RP,nbv,0.d00,gg,nbv)
call dgemm('N','T',nbv,nbv,nbv,1.d00,gg,nbv,T,nbv,0.d00,gg1,nbv)
RP=gg1

  enddo
     call dgetrf(nbv,nbv,RP,nbv,isuppz,info)
     call dgetri(nbv,RP,nbv,isuppz,work,lwork,info)


end subroutine rpropagat

subroutine pott(w,nbv,r,E)
  implicit real(kind=8) (a-h,o-z)
  integer::nbv,l
  real(kind=8),dimension(nbv,nbv),intent(inout) ::w
  real(kind=8)::mu,E,k,r
  mu=1.d0!0.9801045d0!0.9801045d0 !!!!!!!!!!!!µHCl,µCH=0.929931,µCO=6.8606719,µLiH=0.8801221                                  ! Masse reduite
  k=2.d00*mu*E
  l=3


  w(:,:)= (r)*exp(-(r)**2)
  do i=1,nbv
     w(i,i) = -(2.d00*mu*1.d0/(r*r) -  dble(l*(l+1))/((r)**2) + k + dble(i)/100.d0)
  enddo
  return
end subroutine pott
