subroutine ijohnson_d(rmin,rmax,Yh,nbv,n,fl)
!$ use OMP_LIB
        use mod_constantes
        use mod_pot
        use mod_base, only : ltab, jtab, k2tab, lamax, nmax
!        use mod_utildev
  use cublas
  use cusolverDn
  use cudafor

  implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,ir,lam
  integer , intent(in) :: nbv, n
  real(kind=dp),device,allocatable,dimension(:,:)::Qa,Qc,Y,b_d
  real(kind=dp),device,allocatable,dimension(:)::r1,r2
  real(kind=dp),allocatable,dimension(:)::r1h,r2h
  real(kind=dp),allocatable,dimension(:,:)::u
  real(kind=dp),intent(out),dimension(nbv,nbv)::Yh
  real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
  real(kind=dp)::E,k,h,hs3,h2s6,h4,lami
  real(kind=dp), intent(in)::rmin,rmax
  real::result
  real,dimension(2)::tarray
  logical :: iprint
  real(kind=dp) :: BIG=1.d20,one=1.0d0,zero=0.d0

  integer,device,allocatable,dimension(:)::ipiv_d
  integer,device::devinfo_d
  real(kind=dp),allocatable,device,dimension(:)::workspace_d
  type(cusolverDnHandle) :: handle
  TYPE(dim3) :: grid,tblock
!!  cudaDeviceSynchronize()
!!  cudaStreamSynchronize(stream)
!!  cudaStreamQuery(stream)

!!  type (cudaEvent) :: event
!!  ...
!!  istat = cudaEventRecord(event, stream1)
!!  cudaEventSynchronize(event)

  tBlock = dim3(32,4,1)
  grid = dim3(ceiling(real(nbv)/tBlock%x), ceiling(real(nbv)/tBlock%y), 1)

  iprint=.false.
  imp=9
  LIWORK=64*nbv
  h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
  hs3=h/3.d0
  h2s6=-h*h/6.d0
  h4=4.d0/h
  !  call dtime(tarray, result)                                           ! Allocation de la mémoire
  allocate(r1(nbv),r2(nbv))
  allocate(r1h(nbv),r2h(nbv))
  allocate(Qa(nbv,nbv),Qc(nbv,nbv),Y(nbv,nbv))
  allocate(u(nbv,nbv))
  allocate(b_d(nbv,nbv))
  allocate(ipiv_d(nbv))




  istat = cusolverDnCreate(handle)
  istat = cusolverDnDgetrf_bufferSize(handle, nbv, nbv, Y, nbv, Liwork)
  allocate(workspace_d(LiWORK))



  r=rmin
     rm2=1.d0/(r*r)

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
  Y(i,j)=zero
     enddo
     enddo

!$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,nbv
  u(i,j)=zero
     enddo
     enddo
!$OMP END PARALLEL DO
     call pot
  do lam=0,lamax
!$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,nbv
     u(i,j)=u(i,j)+fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!$OMP END PARALLEL DO
  enddo
!$OMP PARALLEL DO PRIVATE (i)
  do i=1,nbv
     u(i,i)=u(i,i)+dfloat(ltab(i)*(ltab(i)+1))*rm2-k2tab(i)
  enddo
!$OMP END PARALLEL DO

Qa=u
!$cuf kernel do(1) <<<*,*>>>
  do i=1,nbv
     Y(i,i)=dsqrt(dabs(Qa(i,i)))
  enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qa(i,i)=zero
     enddo

!     call damat_d<<<grid,tblock>>>(hs3,Qa)
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qa(i,j)=hs3*Qa(i,j)
     enddo
     enddo
!  istat=cudaDeviceSynchronize()
  do ir=1,n,2
     r=r+h                    !*******************  rc
     rm2=1.d0/(r*r)

!$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,nbv
  u(i,j)=zero
     enddo
     enddo
!$OMP END PARALLEL DO
     call pot
  do lam=0,lamax
!$OMP PARALLEL DO PRIVATE (i,j)
     do j=1,nbv
     do i=1,nbv
     u(i,j)=u(i,j)+fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!$OMP END PARALLEL DO
  enddo
!$OMP PARALLEL DO PRIVATE (i)
  do i=1,nbv
     u(i,i)=u(i,i)+dfloat(ltab(i)*(ltab(i)+1))*rm2-k2tab(i)
  enddo
!$OMP END PARALLEL DO

Qc=u

!$cuf kernel do(1) <<<*,*>>>
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

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qc(i,i)=zero
     enddo

!     write(9,*)
!     write(9,*)h2s6,u(1,2),h2s6*u(1,2)
!     write(9,*)h2s6,u(nbv-1,nbv),h2s6*u(nbv-1,nbv)
!     call damat_d<<<grid,tblock>>>(h2s6,Qc)
!u=Qc
!     write(9,*)h2s6,u(1,2)
!     write(9,*)h2s6,u(nbv-1,nbv)

!     call damat_d<<<grid,tblock>>>(h2s6,Qc)
!  istat=cudaDeviceSynchronize()
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qc(i,j)=h2s6*Qc(i,j)
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qc(i,i)=Qc(i,i)+one
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        b_d(i,j)=0.d0
     enddo
     enddo
!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        b_d(i,i)=1.d0
     enddo

  istat = cusolverDnDgetrf(handle, nbv, nbv, Qc, nbv, workspace_d, ipiv_d, devInfo_d)
  istat = cusolverDnDgetrs(handle, CUBLAS_OP_N, nbv, nbv, Qc, nbv, ipiv_d, b_d, nbv, devInfo_d)

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        b_d(i,i)=b_d(i,i)-one
     enddo

!     call damat_d<<<grid,tblock>>>(h4,b_d)
!  istat=cudaDeviceSynchronize()
!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     b_d(i,j)=h4*b_d(i,j)
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qa(i,i)=Qa(i,i)+r1(i)
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Y(i,j)=Y(i,j)+Qa(i,j)
     enddo
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Qc(i,j)=0.d0
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Qc(i,i)=1.d0
     enddo

  istat = cusolverDnDgetrf(handle, nbv, nbv, Y, nbv, workspace_d, ipiv_d, devInfo_d)
  istat = cusolverDnDgetrs(handle, CUBLAS_OP_N, nbv, nbv, Y, nbv, ipiv_d, Qc, nbv, devInfo_d)


     !     call imprime2(iprint,imp,nbv,Y,'Y         ')

!$cuf kernel do(2) <<<*,*>>>
     do i=1,nbv
     do j=1,nbv
        Qc(j,i)=Qc(j,i)*r2(i)
     enddo
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Qc(i,j)=Qc(i,j)*r2(i)
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        b_d(i,i)=b_d(i,i)+r1(i)
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qc(i,j)=b_d(i,j)-Qc(i,j)
     enddo
     enddo

     r=r+h                    !*******************  rb
!$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,nbv
     u(i,j)=zero
     enddo
     enddo
!$OMP END PARALLEL DO
     call pot
     do lam=0,lamax
!$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
     do i=1,nbv
     u(i,j)=u(i,j)+fl(i,j,lam)*vvl(lam)
     enddo
     enddo
!$OMP END PARALLEL DO
     enddo
!$OMP PARALLEL DO PRIVATE (i)
     do i=1,nbv
     u(i,i)=zero
     enddo
!$OMP END PARALLEL DO

     Qa=u

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qa(i,j)=hs3*Qa(i,j)
     enddo
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Qc(i,j)=Qc(i,j)+b_d(i,j)
     enddo
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Y(i,j)=0.d0
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Y(i,i)=1.d0
     enddo

  istat = cusolverDnDgetrf(handle, nbv, nbv, Qc, nbv, workspace_d, ipiv_d, devInfo_d)
  istat = cusolverDnDgetrs(handle, CUBLAS_OP_N, nbv, nbv, Qc, nbv, ipiv_d, Y, nbv, devInfo_d)


!$cuf kernel do(2) <<<*,*>>>
     do i=1,nbv
     do j=1,nbv
        Y(j,i)=Y(j,i)*r2(i)
     enddo
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
        Y(i,j)=Y(i,j)*r2(i)
     enddo
     enddo

!$cuf kernel do(2) <<<*,*>>>
     do j=1,nbv
     do i=1,nbv
     Y(i,j)=Qa(i,j)-Y(i,j)
     enddo
     enddo

!$cuf kernel do(1) <<<*,*>>>
     do i=1,nbv
        Y(i,i)=Y(i,i)+r1(i)
     enddo
!istat=CudaDeviceSynchronize()



  enddo
  Yh=Y

deallocate(ipiv_d)

deallocate(workspace_d)

deallocate(u,r1,r2,Qa,Qc,Y,b_d)
deallocate(r1h,r2h)
istat = devInfo_d
istat = cusolverDnDestroy(handle)


end subroutine ijohnson_d
