subroutine ijohnsonsp(rmin,rmax,Y,nbv,n,nsp)
        use mod_constantes
        use mod_pot
        use mod_base, only : ltab, jtab, k2tab, lamax, nmax
  implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,ir,lam
  integer , intent(in) :: nbv, n, nsp
  real(kind=dp),allocatable,dimension(:,:)::Qa,Qc
  real(kind=dp),allocatable,dimension(:)::QQ
  real(kind=dp),allocatable,dimension(:)::r1,r2
  real(kind=dp),intent(out),dimension(nbv,nbv)::Y
  integer,allocatable,dimension(:)::isuppz
  real(kind=dp),allocatable,dimension(:)::WORK
  real(kind=dp)::E,k,h,hs3,h2s6,h4,lami
  real(kind=dp), intent(in)::rmin,rmax
  real::result
  real,dimension(2)::tarray
  logical :: iprint
  real(kind=dp) :: BIG=1.d20,one=1.0d0,zero=0.d0
  iprint=.false.
  imp=9
  LWORK=64*nbv
  LIWORK=64*nbv
  h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
  hs3=h/3.d0
  h2s6=h*h/6.d0
  h4=4.d0/h
  nbv2=(nbv*(nbv+1))/2
  !  call dtime(tarray, result)                                           ! Allocation de la mémoire
  allocate(r1(nbv),r2(nbv))
  allocate(Qa(nbv,nbv),Qc(nbv,nbv),QQ(nbv2))


  allocate(WORK(LWORK))
  allocate(isuppz(2*nbv))


  r=rmin
     rm2=1.d0/(r*r)

     do j=1,nbv
     do i=1,nbv
  Y(i,j)=zero
     enddo
     enddo

     call pot
     call mkl_cspblas_dcoogemv('N', nbv2, flsp, indexj, indexl, nsp, vvl, QQ)
     do j=1,nbv
     do i=1,j
     ij=(j*(j-1))/2+i
     Qa(i,j)=QQ(ij)
     enddo
     enddo
!          call imprime2(iprint,imp,nbv,Qa,'Qa        ')
  do i=1,nbv
     Qa(i,i)=Qa(i,i)+dble(ltab(i)*(ltab(i)+1))*rm2-k2tab(i)
  enddo

  do i=1,nbv
     Y(i,i)=dsqrt(dabs(Qa(i,i)))
  enddo
     do i=1,nbv
        Qa(i,i)=zero
     enddo

!  w=u
     do j=1,nbv
     do i=1,j
     Qa(i,j)=hs3*Qa(i,j)
     enddo
     enddo

  do ir=1,n,2
     r=r+h                    !*******************  rc
     rm2=1.d0/(r*r)

     call pot
!     vvv(1:lamax+1)=vvl(0:lamax)
     call mkl_cspblas_dcoogemv('N', nbv2, flsp, indexj, indexl, nsp, vvl, QQ)
     do j=1,nbv
     do i=1,j
     ij=(j*(j-1))/2+i
     Qc(i,j)=QQ(ij)
     enddo
     enddo
  do i=1,nbv
     Qc(i,i)=Qc(i,i)+dble(ltab(i)*(ltab(i)+1))*rm2-k2tab(i)
  enddo
       do i=1,nbv
        lami=dsqrt(dabs(Qc(i,i)))
        if (Qc(i,i).GT.0.d00) then
           r1(i)=(lami/dtanh(h*lami))
           r2(i)=(lami/dsinh(h*lami))
        else
           r1(i)=(lami/dtan(h*lami))
           r2(i)=(lami/dsin(h*lami))
        endif
     enddo
     do i=1,nbv
        Qc(i,i)=zero
     enddo

!     Qa=hs3*w
     do j=1,nbv
     do i=1,j
     Qc(i,j)=-h2s6*Qc(i,j)
     enddo
     enddo
!     do j=1,nbv
!     do i=1,j
!     Qc(j,i)=Qc(i,j)
!     enddo
!     enddo
     do i=1,nbv
        Qc(i,i)=Qc(i,i)+one
     enddo
!     call dgetrf(nbv,nbv,Qc,nbv,isuppz,info)
!     call dgetri(nbv,Qc,nbv,isuppz,work,lwork,info)
if (nbv<500)then
     call dsytrf('U',nbv,Qc,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,Qc,nbv,isuppz,work,info)
else
  call matinv(Qc,nbv)
endif
     do i=1,nbv
     do j=i+1,nbv
       Qc(j,i)=Qc(i,j)
     enddo
     enddo

     do i=1,nbv
        Qc(i,i)=Qc(i,i)-one
     enddo
     do j=1,nbv
     do i=1,j
     Qc(i,j)=h4*Qc(i,j)
     enddo
     enddo

     do i=1,nbv
        Qa(i,i)=Qa(i,i)+r1(i)
     enddo
     do j=1,nbv
     do i=1,nbv
     Y(i,j)=Y(i,j)+Qa(i,j)
     enddo
     enddo
!
!

!
!     do j=1,nbv
!     do i=1,j
!     Y(j,i)=Y(i,j)
!     enddo
!     enddo
!     call dgetrf(nbv,nbv,Y,nbv,isuppz,info)
!     call dgetri(nbv,Y,nbv,isuppz,work,lwork,info)
if (nbv<500)then
     call dsytrf('U',nbv,Y,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,Y,nbv,isuppz,work,info)
else
  call matinv(Y,nbv)
endif
     !     call imprime2(iprint,imp,nbv,Y,'Y         ')
!
!

     do i=1,nbv
     do j=1,nbv
        Y(j,i)=Y(j,i)*r2(i)
     enddo
     enddo
     do j=1,nbv
     do i=1,nbv
        Y(i,j)=Y(i,j)*r2(i)
     enddo
     enddo
     do i=1,nbv
        Qc(i,i)=Qc(i,i)+r1(i)
     enddo
     do j=1,nbv
     do i=1,nbv
     Y(i,j)=Qc(i,j)-Y(i,j)
     enddo
     enddo
!
     r=r+h                    !*******************  rb
!
     call pot
!     vvv(1:lamax+1)=vvl(0:lamax)
     call mkl_cspblas_dcoogemv('N', nbv2, flsp, indexj, indexl, nsp, vvl, QQ)
     do j=1,nbv
     do i=1,j
     ij=(j*(j-1))/2+i
     Qa(i,j)=hs3*QQ(ij)
     enddo
     enddo
!
     do i=1,nbv
        Qa(i,i)=zero
     enddo
!
     do j=1,nbv
     do i=1,nbv
     Y(i,j)=Y(i,j)+Qc(i,j)
     enddo
     enddo
     do j=1,nbv
     do i=1,j
     Y(j,i)=Y(i,j)
     enddo
     enddo
if (nbv<500)then
     call dsytrf('U',nbv,Y,nbv,isuppz,work,lwork,info)
     call dsytri('U',nbv,Y,nbv,isuppz,work,info)
else
  call matinv(Y,nbv)
endif

     do i=1,nbv
     do j=1,nbv
        Y(j,i)=Y(j,i)*r2(i)
     enddo
     enddo
     do j=1,nbv
     do i=1,nbv
        Y(i,j)=Y(i,j)*r2(i)
     enddo
     enddo
     do j=1,nbv
     do i=1,nbv
     Y(i,j)=Qa(i,j)-Y(i,j)
     enddo
     enddo
     do i=1,nbv
        Y(i,i)=Y(i,i)+r1(i)
     enddo

     do j=1,nbv
     do i=j,nbv
     Y(i,j)=Y(j,i)
     enddo
     enddo


  enddo
!          call imprime2(iprint,imp,nbv,Y,'Y         ')
deallocate(r1,r2,Qa,Qc,QQ,work,isuppz)

end subroutine ijohnsonsp
