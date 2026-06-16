subroutine tmat(RP,nbv,nopen,ltab,ktab,rmax,Tr,Ti,iprint,imp)
  !$ use OMP_LIB
  use mod_constantes, only : pi,dp
  implicit real(kind=dp) (a-h,o-z)
  !integer, parameter :: dp = 8
  !integer, parameter :: dp = selected_real_kind(15, 307)
  integer, intent(in) :: nbv,nopen,  iprint, imp
  integer, intent(in), dimension(nbv) :: ltab
  real(kind=dp), intent(in), dimension(nbv) :: ktab
  real(kind=dp), intent(in), dimension(nbv,nbv) ::RP
  real(kind=dp), intent(in) :: rmax
  real(kind=dp), intent(out), dimension(nopen,nopen) :: Tr,Ti
  real(kind=dp), allocatable,  dimension(:,:) :: Jn,Nn,Jnp,Nnp,K,RPo,Ko
  real(kind=dp), allocatable, dimension(:) :: work
  real(kind=dp) :: z, k12, fac, a1, a2, rl
  integer ::info,lwork,irest,i,j
  integer, allocatable, dimension(:) :: ipiv


  allocate(ipiv(nbv))
  allocate(Jn(nbv,nbv),Nn(nbv,nbv),Jnp(nbv,nbv),Nnp(nbv,nbv))
  allocate(K(nbv,nbv),Ko(nopen,nopen),RPo(nopen,nopen))

  RPo=RP(1:nopen,1:nopen)
  !call imprime2(iprint,imp,nopen,RPo,'RPopen    ')
  lwork=64*nbv
  allocate(work(lwork))
  Jn=0.d0
  Jnp=0.d0
  Nn=0.d0
  Nnp=0.d0
  K=0.d0
  !$OMP PARALLEL
  !$OMP DO PRIVATE (i,k12,z,sqk2,AJn,ANn,fac)
  do i=1,nbv
     k12=dsqrt(dabs(ktab(i)))
     if(ktab(i)<0.d0)k12=-k12
     z=rmax*k12
     !if(ktab(i)<0.d0) then
     !call besselm(ltab(i),z,Jn(i,i),Jnp(i,i),Nn(i,i),Nnp(i,i))
     !else
     call bessel(ltab(i),z,Jn(i,i),Jnp(i,i),Nn(i,i),Nnp(i,i))
     if(ktab(i)>0.d0)then
        k12=dabs(k12)
        z=dabs(z)
        sqk12=1.d0/dsqrt(k12)
        AJn=Jn(i,i)
        Jn(i,i)=sqk12*z*AJn
        Jnp(i,i)=sqk12*k12*(AJn+z*Jnp(i,i))
        ANn=Nn(i,i)
        Nn(i,i)=sqk12*z*ANn
        Nnp(i,i)=sqk12*k12*(ANn+z*Nnp(i,i))
     else
        k12=dabs(k12)
        z=dabs(z)
        fac=dsqrt(pi/(2.d0*z))
        AJn=Jn(i,i)
        Jn(i,i)=Jn(i,i)*fac
        !Jnp(i,i)=fac*(Jnp(i,i)-AJn*k12/(2.d0*z))
        Jnp(i,i)=fac*k12*(Jnp(i,i)-AJn/(2.d0*z))
        ANn=Nn(i,i)
        Nn(i,i)=Nn(i,i)*fac
        !Nnp(i,i)=fac*(Nnp(i,i)-ANn*k12/(2.d0*z))
        Nnp(i,i)=fac*k12*(Nnp(i,i)-ANn/(2.d0*z))
        !Jn(i,i)=0.d0
        !Nn(i,i)=0.d0
        !Jnp(i,i)=5.d9
        !Nnp(i,i)=5.d9
     endif
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  !call imprime2(iprint,imp,nopen,Jn,'Jn        ')
  !call imprime2(iprint,imp,nopen,Jnp,'Jnp       ')
  !call imprime2(iprint,imp,nopen,Nn,'Nn        ')
  !call imprime2(iprint,imp,nopen,Nnp,'Nnp       ')
  !
  !                  matrice K
  !
  call dgemm('N','N',nbv,nbv,nbv,-1.d00,RP,nbv,Jn,nbv,1.d00,Jnp,nbv)
  !call imprime2(iprint,imp,nopen,Jnp,'Jn-RP*Jnp ')
  call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,Nn,nbv,-1.d00,Nnp,nbv)
  !ggr=Nnp
  call dgetrf(nbv,nbv,Nnp,nbv,ipiv,info)
  call dgetri(nbv,Nnp,nbv,ipiv,work,lwork,info)
  !call imprime2(iprint,imp,nopen,Nnp,'Num(-1)   ')
  call dgemm('N','N',nbv,nbv,nbv,1.d00,Nnp,nbv,Jnp,nbv,0.d00,K,nbv)
  !call imprime2(iprint,imp,nopen,K,'K         ')
  Ko=K(1:nopen,1:nopen)
  !call imprime2(iprint,imp,nopen,Ko,'Ko        ')

  do i=1,nopen
     do j=1,i-1
        Ko(i,j)=(Ko(i,j)+Ko(j,i))/2.d0
        Ko(j,i)=Ko(i,j)
     enddo
  enddo
  !call imprime2(iprint,imp,nopen,Ko,'Ko        ')
  !
  !--------------------------------------------------
  !                  matrice T = Tr + i Ti
  !
  !LW   relation between matrices (Newton, eq. 7.57/58/59:
  !  pi*T=-K*(1-iK)^(-1)
  !  S = 1 - 2*pi*i*T
  !  S=(1+iK)*(1-i*K)^(-1)
  !--------------------------------------------------

  call dgemm('N','N',nopen,nopen,nopen,1.d00,Ko,nopen,Ko,nopen,0.d00,Tr,nopen)
  do i=1,nopen
     Tr(i,i)=1.d0+Tr(i,i)
  enddo
  !call imprime2(iprint,imp,nopen,Tr,'I+K2      ')
  call dgetrf(nopen,nopen,Tr,nopen,ipiv,info)
  call dgetri(nopen,Tr,nopen,ipiv,work,lwork,info)
  !call imprime2(iprint,imp,nopen,Tr,'T-        ')
  call dgemm('N','N',nopen,nopen,nopen,-2.d00,Ko,nopen,Tr,nopen,0.d00,Ti,nopen)
  call dgemm('N','N',nopen,nopen,nopen,-1.d00,Ti,nopen,Ko,nopen,0.d00,Tr,nopen)
  !call imprime2(iprint,imp,nopen,Tr,'Tr        ')
  !call imprime2(iprint,imp,nopen,Ti,'Ti        ')
  !
  deallocate(ipiv)
  deallocate(Jn,Nn,Jnp,Nnp)
  deallocate(K,Ko,RPo)
end subroutine tmat
