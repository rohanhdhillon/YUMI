subroutine johnson(rmin,rmax,fl,Y,nbv,n)
        use mod_constantes
        use mod_pot
        use mod_base, only : ltab, jtab, k2tab, lamax, nmax
!$ use OMP_LIB
  implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,ir,lam,i,j
  integer , intent(in) :: nbv, n
  real(kind=dp),allocatable,dimension(:,:)::w,u,T1
  real(kind=dp),intent(out),dimension(nbv,nbv)::Y
  real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
  integer,allocatable,dimension(:)::isuppz
  real(kind=dp),allocatable,dimension(:)::WORK
  real(kind=dp)::E,k,h,h2s3,h2s6,m4h2s3,m2h2s3
  real(kind=dp), intent(in)::rmin,rmax
  logical :: iprint
  real(kind=dp) :: BIG=1.d20
  iprint=.false.
  !LWORK=-1
  !LIWORK=-1
  LWORK=64*nbv
  LIWORK=64*nbv
  h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
  h2s3=h*h/3.d0
  h2s6=h2s3/2.d0
  m4h2s3=-4.d0*h2s3
  m2h2s3=-2.d0*h2s3
  !                                                            ! Allocation de la mémoire
  allocate(w(nbv,nbv))
  allocate(T1(nbv,nbv))
  allocate(u(nbv,nbv))


  allocate(WORK(LWORK))
  allocate(isuppz(nbv))


  r=rmin
  rm2=1.d0/(r*r)

  Y(:,:)=0.0d00
  do i=1,nbv
     Y(i,i)=BIG
  enddo


     u(:,:)=0.d0
     call pot
     do lam=0,lamax
     do j=1,nbv
     do i=j,nbv
        u(i,j)=u(i,j)-fl(i,j,lam)*vvl(lam)
     enddo
     enddo
     enddo
     do i=1,nbv
        u(i,i)=u(i,i)-dfloat(ltab(i)*(ltab(i)+1))*rm2+k2tab(i)
     enddo
     do i=1,nbv
     do j=i+1,nbv
        u(i,j)=u(j,i)
     enddo
     enddo

     T1=Y
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,Y,nbv,-h2s3,u,nbv)
     Y=u
     r=r+h
  rm2=1.d0/(r*r)

  do ir=1,n-3,2
     w(:,:)=0.d0
     call pot
     do lam=0,lamax
     do j=1,nbv
     do i=j,nbv
        w(i,j)=w(i,j)-fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!        w(:,:)=w(:,:)-fl(:,:,lam)*vvl(lam)
     enddo
     do i=1,nbv
        w(i,i)=w(i,i)-dfloat(ltab(i)*(ltab(i)+1))*rm2+k2tab(i)
     enddo
     do i=1,nbv
     do j=i+1,nbv
        w(i,j)=w(j,i)
     enddo
     enddo

     T1=h2s6*w
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,w,nbv,0.d0,u,nbv)
     T1=Y
!  call imprime2(iprint,2,nbv,Y,'Y         ')
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,Y,nbv,m4h2s3,u,nbv)
     Y=u
     r=r+h
  rm2=1.d0/(r*r)

     u(:,:)=0.d0
     call pot
     do lam=0,lamax
     do j=1,nbv
     do i=j,nbv
        u(i,j)=u(i,j)-fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!        u(:,:)=u(:,:)-fl(:,:,lam)*vvl(lam)
     enddo
     do i=1,nbv
        u(i,i)=u(i,i)-dfloat(ltab(i)*(ltab(i)+1))/(r*r)+k2tab(i)
     enddo
     do i=1,nbv
     do j=i+1,nbv
        u(i,j)=u(j,i)
     enddo
     enddo

     T1=Y
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,Y,nbv,m2h2s3,u,nbv)
     Y=u
     r=r+h


  enddo


     w(:,:)=0.d0
     call pot
     do lam=0,lamax
     do j=1,nbv
     do i=j,nbv
        w(i,j)=w(i,j)-fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!        w(:,:)=w(:,:)-fl(:,:,lam)*vvl(lam)
     enddo
     do i=1,nbv
        w(i,i)=w(i,i)-dfloat(ltab(i)*(ltab(i)+1))/(r*r)+k2tab(i)
     enddo
     do i=1,nbv
     do j=i+1,nbv
        w(i,j)=w(j,i)
     enddo
     enddo

     T1=h2s6*w
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,w,nbv,0.d0,u,nbv)
     T1=Y
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,Y,nbv,m4h2s3,u,nbv)
     Y=u
     r=r+h
  rm2=1.d0/(r*r)

     u(:,:)=0.d0
     call pot
     do lam=0,lamax
     do j=1,nbv
     do i=j,nbv
        u(i,j)=u(i,j)-fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!        u(:,:)=u(:,:)-fl(:,:,lam)*vvl(lam)
     enddo
     do i=1,nbv
        u(i,i)=u(i,i)-dfloat(ltab(i)*(ltab(i)+1))*rm2+k2tab(i)
     enddo
     do i=1,nbv
     do j=i+1,nbv
        u(i,j)=u(j,i)
     enddo
     enddo

     T1=Y
     do i=1,nbv
         T1(i,i)=1.d0+T1(i,i)
     enddo
!     call dgetrf(nbv,nbv,T1,nbv,isuppz,info)
!     call dgetri(nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytrf('U',nbv,T1,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,T1,nbv,isuppz,work,info)
     do i=1,nbv
     do j=i+1,nbv
       T1(j,i)=T1(i,j)
     enddo
     enddo

     call dgemm('N','N',nbv,nbv,nbv,1.d00,T1,nbv,Y,nbv,-h2s3,u,nbv)
     Y=u/h
!write(*,*) 'rmax=',r

end subroutine johnson











