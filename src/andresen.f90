subroutine andresen(rmin,rmax,ri,fi,fs,hh,fl,nmax,Y,nbv,lamax,jtab,ltab,ktab,n)
  implicit real(kind=8) (a-h,o-z)
integer, parameter :: dp = selected_real_kind(15, 307)
  integer::INFO,LWORK,LIWORK,ir,lam
  real(kind=dp),allocatable,dimension(:,:)::w,u,T
  real(kind=dp),allocatable,dimension(:,:)::Y1,Y1p
  real(kind=dp),intent(out),dimension(nbv,nbv)::Y
  real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
  integer, intent(in), dimension(nbv) :: jtab, ltab
  real(kind=dp), intent(in), dimension(nbv) :: ktab
  integer,allocatable,dimension(:)::isuppz,iwork
  real(kind=dp),allocatable,dimension(:)::WORK,C11,C12,C21,C22,phi,phip,phi1
  real(kind=dp),allocatable,dimension(:)::Jn,Jnp,Nn,Nnp,Jn1,Jn1p,Nn1,Nn1p,lambda,q
  real(kind=dp)::E,k,r,h,h2s3,h2s6,m4h2s3,m2h2s3,pi,vl
  integer , intent(in) :: nbv, lamax,n
  real(kind=dp),intent(in), dimension(0:nmax) :: ri,hh
  real(kind=dp),intent(in), dimension(0:nmax,0:lamax) :: fi,fs
  real(kind=dp), intent(in)::rmin,rmax
  real(kind=dp)::vvl
  real::result
  real,dimension(2)::tarray
  logical :: iprint
  real(kind=dp) :: BIG=1.d20
  iprint=.false.
  pi=4.d0*datan(1.d0)
  lda=nbv
  ldz=nbv
  !LWORK=-1
  !LIWORK=-1
  LWORK=64*nbv
  LIWORK=64*nbv
  !n=400
  h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
  hs2=h/2.d0
  call dtime(tarray, result)                                           ! Allocation de la mémoire
  allocate(w(nbv,nbv))
  allocate(T(nbv,nbv),Y1(nbv,nbv),Y1p(nbv,nbv))
  allocate(u(nbv,nbv))
  allocate(lambda(nbv),q(nbv))
  allocate(C11(nbv),C12(nbv),C21(nbv),C22(nbv),phi(nbv),phip(nbv),phi1(nbv))
  allocate(Jn(nbv),Jnp(nbv),Nn(nbv),Nnp(nbv),Jn1(nbv),Jn1p(nbv),Nn1(nbv),Nn1p(nbv))


  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  allocate(isuppz(2*nbv))

  r=rmin
     rm=r+hs2
     w(:,:)=0.d0
     do lam=0,lamax
        call pot(ri,fi(:,lam),fs(:,lam),hh,nmax,rm,vvl)
        w(:,:)=w(:,:)+fl(:,:,lam)*vvl
     enddo

  call DSYEVR('V','A','L',nbv,w,lda,vl,vu,il,iu,abstol,M,lambda, &
       T,ldz,isuppz,work,LWORK,iwork,LiWORK,info) 

phip=h*(-dsqrt(dabs(lambda))+h*lambda/3.d0)
phi=1.d0


  do ir=1,n-1
     rm=r+hs2
     w(:,:)=0.d0
     do lam=0,lamax
        call pot(ri,fi(:,lam),fs(:,lam),hh,nmax,rm,vvl)
        w(:,:)=w(:,:)+fl(:,:,lam)*vvl
     enddo

  call DSYEVR('V','A','L',nbv,w,lda,vl,vu,il,iu,abstol,M,lambda, &
       T,ldz,isuppz,work,LWORK,iwork,LiWORK,info) 

!write(9,*)'Lambda=', lambda

do i=1,nbv
if((ktab(i)-lambda(i))*r>=0.d0)then
k12=dsqrt(ktab(i)-lambda(i))
z=k12*r1
!z=ktab(i)-lambda(i)*r
call bessel(ltab(i),z,Jn(i),Jnp(i),Nn(i),Nnp(i))
else
k12=dsqrt(-ktab(i)+lambda(i))
z=k12*r1
!z=+(ktab(i)-lambda(i)*r)
Jn(i)=1.d0
Nn(i)=1.d0
a1=1.d0
do j=1,ltab(i)
a2=a1+dfloat(j+j-1)/z
a1=1.d0/a2
enddo
rl=dfloat(ltab(i))
Nnp(i)=(1.d0/z-(rl*a1+(rl+1.d0)*a2)/(rl+rl+1.d0))*z
Jnp(i)=Nnp(i)
endif
enddo

do i=1,nbv
r1=r+h
if((ktab(i)-lambda(i))*r1>=0.d0)then
k12=dsqrt(ktab(i)-lambda(i))
!z=ktab(i)-lambda(i)*r1
z=k12*r1
call bessel(ltab(i),z,Jn1(i),Jn1p(i),Nn1(i),Nn1p(i))
else
k12=dsqrt(-ktab(i)+lambda(i))
!z=+(ktab(i)-lambda(i)*r1)
z=k12*r1
Jn1(i)=1.d0
Nn1(i)=1.d0
a1=1.d0
do j=1,ltab(i)
a2=a1+dfloat(j+j-1)/z
a1=1.d0/a2
enddo
rl=dfloat(ltab(i))
Nn1p(i)=(1.d0/z-(rl*a1+(rl+1.d0)*a2)/(rl+rl+1.d0))*z
Jn1p(i)=Nnp(i)
endif
enddo

q=dsqrt(dabs(ktab-lambda))
C11=Jn1*Nnp-Nn1*Jnp
C12=(Nn1*Jn-Jn1*Nn)/q
C21=(Jn1p*Nnp-Nn1p*Jnp)*q
C22=Nn1p*Jn-Jn1p*Nn

phi1=C11*phi+C12*phip
phip=C21*phi+C22*phip
phi=phi1


r=r+h

  enddo

phi1=phip/phi

do i=1,nbv
do j=1,nbv
   Y(i,j)=T(j,i)*phi1(i)
!   Y1(i,j)=T(j,i)*phi(j)
!   Y1p(i,j)=T(j,i)*phip(j)
enddo
enddo
!call dgemm('N','N',nbv,nbv,nbv,1.d00,Y1p,nbv,Y1,nbv,0.d00,Y,nbv)

end subroutine andresen
