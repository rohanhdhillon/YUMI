subroutine numerov(rmin,rmax,ri,fi,fs,hh,fl,nmax,Y,nbv,lamax,jtab,ltab,ktab,n,leftright,h,noe)
  implicit real(kind=8) (a-h,o-z)
  integer, parameter :: dp = 8
  !integer, parameter :: dp = selected_real_kind(15, 307)
  integer::INFO,LWORK,LIWORK,ir,lam
  integer , intent(in) :: nbv, lamax,n
  integer , intent(inout) :: noe
  real(kind=dp),dimension(nbv,nbv)::w,u,T1,Y2
  real(kind=dp),intent(out),dimension(nbv,nbv)::Y
  real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
  integer, intent(in), dimension(nbv) :: jtab, ltab
  real(kind=dp), intent(in), dimension(nbv) :: ktab
  real(kind=dp), dimension(nbv) :: vp
  integer,allocatable,dimension(:)::isuppz
  real(kind=dp),allocatable,dimension(:)::WORK
  real(kind=dp)::E,k,r,h2s3,h2s6,m4h2s3,m2h2s3,pi,vvl,h2s12
  real(kind=dp),intent(in), dimension(0:nmax) :: ri,hh
  real(kind=dp),intent(in), dimension(0:nmax,0:lamax) :: fi,fs
  real(kind=dp), intent(in)::rmin,rmax,leftright,h
  real::result
  real,dimension(2)::tarray
  logical :: iprint
  real(kind=dp) :: BIG=1.d40
  iprint=.false.
  pi=4.d0*datan(1.d0)
  !LWORK=-1
  !LIWORK=-1
  LWORK=64*nbv
  LIWORK=64*nbv
  !n=400
!  h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
  h2s3=h*h/3.d0
  h2s12=h*h/12.d0
  h2s6=h2s3/2.d0
  m4h2s3=-4.d0*h2s3
  m2h2s3=-2.d0*h2s3
!  call dtime(tarray, result)                                           ! Allocation de la mémoire
!  allocate(w(nbv,nbv))
!  allocate(T1(nbv,nbv))
!  allocate(Y2(nbv,nbv))
!  allocate(u(nbv,nbv))


  allocate(WORK(LWORK))
  allocate(isuppz(2*nbv))

  if( leftright.GT.0.d0) then
     r=rmin
  else
     r=rmax
  endif

  Y2(:,:)=0.0d00
!  noe=0


  do ir=1,n
     u(:,:)=0.d0
     do lam=0,lamax
        call pot(ri,fi(:,lam),fs(:,lam),hh,nmax,r,vvl)
        u(:,:)=u(:,:)-fl(:,:,lam)*vvl
     enddo
     do i=1,nbv
        u(i,i)=u(i,i)-dfloat(ltab(i)*(ltab(i)+1))/(r*r)+ktab(i)
     enddo

!     write(2,*) 'r=',r
!     call imprime2(iprint,2,nbv,u,'u         ')
     
     T1=h2s12*u
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call imprime2(iprint,2,nbv,T1,'T1        ')
     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     T1=12.d0*T1
     do i=1,nbv
         T1(i,i)=T1(i,i)-10.d0
     enddo
     Y=T1-Y2
     Y2=Y
!     call imprime2(iprint,2,nbv,Y,'Y         ')
     call dgetrf(nbv,nbv,Y2,nbv,isuppz,info)
     prodY=1.0d0
     do i=1,nbv
        if( Y2(i,i).LT.0.d0) noe=noe+1
        prodY=prodY*Y2(i,i)
!        write(*,*) i, Y2(i,i)
     enddo
!        if( prodY.LT.0.d0) noe=noe+1
     call dgetri(nbv,Y2,nbv,isuppz,work,lwork,info)

     r=r+leftright*h
!        write(*,*) ir,r,u(1,1)
  enddo

  if( leftright.LT.0.d0) then
     Y=Y2  
  endif


!write(*,*) 'noe=', noe
!pause
end subroutine numerov
