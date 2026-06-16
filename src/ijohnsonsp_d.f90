subroutine ijohnsonsp_d(rmin,rmax,Y,nbv,n,nsp)
!$ use OMP_LIB
        use mod_constantes
        use mod_pot
        use mod_base, only : ltab, jtab, k2tab, lamax, nmax
  use cudafor
  use cublas
  use cusolverDn
  implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,ir,lam
  integer , intent(in) :: nbv, n, nsp
  real(kind=dp),device,allocatable,dimension(:,:)::Qa_d,Qc_d,Y_d,Qb_d
  real(kind=dp),device,allocatable,dimension(:)::r1_d,r2_d
  real(kind=dp),device,allocatable,dimension(:)::k2tab_d,ltab_d
  real(kind=dp),allocatable,dimension(:,:)::Qa,Qc,Qb
  real(kind=dp),allocatable,dimension(:)::QQ
  real(kind=dp),allocatable,dimension(:)::r1,r2
  real(kind=dp),intent(out),dimension(nbv,nbv)::Y
  integer,allocatable,dimension(:)::isuppz
  real(kind=dp),allocatable,dimension(:)::WORK
  real(kind=dp)::E,k,h,hs3,h2s6,h4,lami
!  real(kind=dp)::E,k,h,hs3,h2s6,h4,lami
!  real(kind=dp),device::h,hs3,h2s6,h4
  real(kind=dp), intent(in)::rmin,rmax
  real::result
  real,dimension(2)::tarray
  logical :: iprint
!  real(kind=dp) :: BIG=1.d20,one=1.0d0,zero=0.d0
  real(kind=dp) :: BIG=1.d20
  real(kind=dp),device :: one,zero



  integer,device,allocatable,dimension(:)::ipiv_d
  integer,device::devinfo_d
  real(kind=dp),allocatable,device,dimension(:)::workspace_d
  type(cusolverDnHandle) :: handle
  TYPE(dim3) :: grid,tblock



  tBlock = dim3(32,4,1)
  grid = dim3(ceiling(real(nbv)/tBlock%x), ceiling(real(nbv)/tBlock%y), 1)

  one=1.d0
  zero=0.d0

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
  allocate(r1_d(nbv),r2_d(nbv))
  allocate(k2tab_d(nbv),ltab_d(nbv))
  allocate(Qa_d(nbv,nbv),Qc_d(nbv,nbv),Y_d(nbv,nbv),Qb_d(nbv,nbv))


  allocate(WORK(LWORK))
  allocate(isuppz(2*nbv))
   
  allocate(ipiv_d(nbv))

  istat = cusolverDnCreate(handle)
  istat = cusolverDnDgetrf_bufferSize(handle, nbv, nbv, Y_d, nbv, Liwork)
  allocate(workspace_d(LiWORK))


  k2tab_d=k2tab
  ltab_d=ltab

  r=rmin
     rm2=1.d0/(r*r)

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
  Y_d(i,j)=zero
     enddo
     enddo

     call pot
     call mkl_cspblas_dcoogemv('N', nbv2, flsp, indexj, indexl, nsp, vvl, QQ)
!$OMP PARALLEL
!$OMP DO PRIVATE (i,j,ij)
     do j=1,nbv
     do i=1,j
     ij=(j*(j-1))/2+i
     Qa(i,j)=QQ(ij)
     enddo
     enddo
!$OMP END DO
!$OMP DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,j
     Qa(j,i)=Qa(i,j)
     enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
!          call imprime2(iprint,imp,nbv,Qa,'Qa        ')
Qa_d=Qa
!$cuf kernel do(1) <<<*,*>>>
  do i=1,nbv
     Qa_d(i,i)=Qa_d(i,i)+dble(ltab_d(i)*(ltab_d(i)+1))*rm2-k2tab_d(i)
  enddo

!$cuf kernel do(1) <<<*,*>>>
  do i=1,nbv
     Y_d(i,i)=dsqrt(dabs(Qa_d(i,i)))
  enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qa_d(i,i)=zero
     enddo

!  w=u
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qa_d(i,j)=hs3*Qa_d(i,j)
     enddo
     enddo
!  Qa=Qa_d
!  Y=Y_d

  do ir=1,n,2
     r=r+h                    !*******************  rc
     rm2=1.d0/(r*r)
 ! Qa=Qa_d
 ! Y=Y_d

     call pot
!     vvv(1:lamax+1)=vvl(0:lamax)
     call mkl_cspblas_dcoogemv('N', nbv2, flsp, indexj, indexl, nsp, vvl, QQ)
!$OMP PARALLEL
!$OMP DO PRIVATE (i,j,ij)
     do j=1,nbv
     do i=1,j
     ij=(j*(j-1))/2+i
     Qc(i,j)=QQ(ij)
     enddo
     enddo
!$OMP END DO
!$OMP DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,j
     Qc(j,i)=Qc(i,j)
     enddo
     enddo
!$OMP END DO
!$OMP END PARALLEL
Qc_d=Qc
!$cuf kernel do(1) <<<*,*>>>
  do i=1,nbv
     Qc_d(i,i)=Qc_d(i,i)+dble(ltab_d(i)*(ltab_d(i)+1))*rm2-k2tab_d(i)
  enddo
!$cuf kernel do(1) <<<*,*>>>
       do i=1,nbv
        lami=dsqrt(dabs(Qc_d(i,i)))
        if (Qc_d(i,i).GT.0.d00) then
           r1_d(i)=(lami/dtanh(h*lami))
           r2_d(i)=(lami/dsinh(h*lami))
        else
           r1_d(i)=(lami/dtan(h*lami))
           r2_d(i)=(lami/dsin(h*lami))
        endif
     enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qc_d(i,i)=zero
     enddo

!     Qa=hs3*w
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qc_d(i,j)=-h2s6*Qc_d(i,j)
     enddo
     enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qc_d(i,i)=Qc_d(i,i)+one
     enddo
!Qc=Qc_d
!     call dgetrf(nbv,nbv,Qc,nbv,isuppz,info)
!     call dgetri(nbv,Qc,nbv,isuppz,work,lwork,info)


!if (nbv<500)then
!     call dsytrf('U',nbv,Qc,nbv,isuppz,work,lwork,info)
!     call dsytri('U',nbv,Qc,nbv,isuppz,work,info)
!else
!  call matinv(Qc,nbv)
!endif

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Qb_d(i,j)=zero
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qb_d(i,i)=one
     enddo


  istat = cusolverDnDgetrf(handle, nbv, nbv, Qc_d, nbv, workspace_d, ipiv_d, devInfo_d)
  istat = cusolverDnDgetrs(handle, CUBLAS_OP_N, nbv, nbv, Qc_d, nbv, ipiv_d, Qb_d, nbv, devInfo_d)



!     do i=1,nbv
!!!$OMP PARALLEL DO PRIVATE (j )
!     do j=i+1,nbv
!       Qc(j,i)=Qc(i,j)
!     enddo
!!!$OMP END PARALLEL DO
!     enddo

!!!!!$OMP PARALLEL DO PRIVATE (i , j )
!     do j=1,nbv
!     do i=1,j
!     Qc(j,i)=Qc(i,j)
!     enddo
!     enddo
!!!!!$OMP END PARALLEL DO
!Qc_d=Qb_d
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qb_d(i,i)=Qb_d(i,i)-one
     enddo
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qb_d(i,j)=h4*Qb_d(i,j)
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qa_d(i,i)=Qa_d(i,i)+r1_d(i)
     enddo
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Y_d(i,j)=Y_d(i,j)+Qa_d(i,j)
     enddo
     enddo
!
!

!
!Y=Y_d
!     call dgetrf(nbv,nbv,Y,nbv,isuppz,info)
!     call dgetri(nbv,Y,nbv,isuppz,work,lwork,info)
!if (nbv<500)then
!     call dsytrf('U',nbv,Y,nbv,isuppz,work,lwork,info)
!     call dsytri('U',nbv,Y,nbv,isuppz,work,info)
!else
!  call matinv(Y,nbv)
!endif
     !     call imprime2(iprint,imp,nbv,Y,'Y         ')
!
!
!!$OMP PARALLEL DO PRIVATE (i , j )
!     do j=1,nbv
!     do i=1,j
!     Y(j,i)=Y(i,j)
!     enddo
!     enddo
!!$OMP END PARALLEL DO
!Y_d=Y


!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Qc_d(i,j)=zero
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qc_d(i,i)=one
     enddo


  istat = cusolverDnDgetrf(handle, nbv, nbv, Y_d, nbv, workspace_d, ipiv_d, devInfo_d)
  istat = cusolverDnDgetrs(handle, CUBLAS_OP_N, nbv, nbv, Y_d, nbv, ipiv_d, Qc_d, nbv, devInfo_d)

 
!Y_d=Qc_d
!$cuf kernel do(2) <<<*,*>>>
     do i=1,nbv
     do j=1,nbv
        Qc_d(j,i)=Qc_d(j,i)*r2_d(i)
     enddo
     enddo
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Qc_d(i,j)=Qc_d(i,j)*r2_d(i)
     enddo
     enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qb_d(i,i)=Qb_d(i,i)+r1_d(i)
     enddo
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qc_d(i,j)=Qb_d(i,j)-Qc_d(i,j)
     enddo
     enddo
!
     r=r+h                    !*******************  rb
!
     call pot
!     vvv(1:lamax+1)=vvl(0:lamax)
     call mkl_cspblas_dcoogemv('N', nbv2, flsp, indexj, indexl, nsp, vvl, QQ)
!$OMP PARALLEL DO PRIVATE (i,j,ij)
     do j=1,nbv
     do i=1,j
     ij=(j*(j-1))/2+i
     Qa(i,j)=hs3*QQ(ij)
     enddo
     enddo
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,j
     Qa(j,i)=Qa(i,j)
     enddo
     enddo
!$OMP END PARALLEL DO
Qa_d=Qa
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qa_d(i,i)=zero
     enddo
!
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qc_d(i,j)=Qc_d(i,j)+Qb_d(i,j)
     enddo
     enddo
!Y=Y_d
!if (nbv<500)then
!     call dsytrf('U',nbv,Y,nbv,isuppz,work,lwork,info)
!     call dsytri('U',nbv,Y,nbv,isuppz,work,info)
!else
!  call matinv(Y,nbv)
!endif
!!$OMP PARALLEL DO PRIVATE (i , j )
!     do j=1,nbv
!     do i=1,j
!     Y(j,i)=Y(i,j)
!     enddo
!     enddo
!!$OMP END PARALLEL DO
!Y_d=Y




!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Y_d(i,j)=zero
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Y_d(i,i)=one
     enddo


  istat = cusolverDnDgetrf(handle, nbv, nbv, Qc_d, nbv, workspace_d, ipiv_d, devInfo_d)
  istat = cusolverDnDgetrs(handle, CUBLAS_OP_N, nbv, nbv, Qc_d, nbv, ipiv_d, Y_d, nbv, devInfo_d)


!Y_d=Qb_d

!$cuf kernel do(2) <<<*,*>>>
     do i=1,nbv
     do j=1,nbv
        Y_d(j,i)=Y_d(j,i)*r2_d(i)
     enddo
     enddo
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Y_d(i,j)=Y_d(i,j)*r2_d(i)
     enddo
     enddo
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Y_d(i,j)=Qa_d(i,j)-Y_d(i,j)
     enddo
     enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Y_d(i,i)=Y_d(i,i)+r1_d(i)
     enddo



  enddo
Y=Y_d
!          call imprime2(iprint,imp,nbv,Y,'Y         ')
deallocate(r1,r2,Qa,Qc,QQ,work,isuppz,Qa_d,Qc_d,r1_d,r2_d,ltab_d,k2tab_d,Qb_d)

end subroutine ijohnsonsp_d
