subroutine V03(jtot,nlev,fl)
!  Lecture du potentiel
  !  Version 1.5
  !NJ   26/4/22  + LW
  !
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
!$ use OMP_LIB
  implicit none
  integer :: i, j, k,  ndim, lam, irest
  integer :: lmin,lmax
  integer, intent(in) :: jtot,nlev
  real(kind=dp) :: fac,som, f3j0,faci,facj,facl
  real(kind=dp), intent(out), dimension(nlev,nlev,0:lamax) :: fl
  logical :: printw
  real(kind=qp) :: w3js, w6js, w9js, som1, ninej,sixj,threej0,som2,som3,som4
  real(kind=qp) :: aa1,aa2,pi32,eps
  integer :: j1i,j2i,j12i,llli,j1j,j2j,j12j,lllj,ll2,il,ll1,ll,ipar
  ! write(9,*) fact
  printw=.false.

  aa1=pi**(1.5d0)
  pi32=(4.d0*pi)**(-1.5d0)

  aa2=5.56832799683171d0
  eps=1.d-10
  
  
  !$OMP PARALLEL DO PRIVATE(i,j,il)
  do j=1,nlev
     do il=0,lamax
        do i=1,nlev
            fl(i,j,il)=0.d0
        enddo
     enddo
  enddo
  !OMP END PARALLEL DO

  
  
  
  !$OMP PARALLEL DO PRIVATE(il,som1,som2,som3,som4,som,facl,ll,ll1,ll2,j1i,j2i,j12i,j1j,j2j,j12j,llli,lllj,facj,faci,ipar)
  do i=1,nlev
      j1i=j1tab(i)
      j2i=j2tab(i)
      j12i=j12tab(i)
      llli=ltab(i)
      faci=dsqrt(dfloat((2*j1i+1)*(2*j2i+1)*(2*j12i+1)*(2*llli+1)))*pi32
      do j=1,nlev
         j1j=j1tab(j)
         j2j=j2tab(j)
         j12j=j12tab(j)
         lllj=ltab(j)
         facj=dsqrt(dfloat((2*j1j+1)*(2*j2j+1)*(2*j12j+1)*(2*lllj+1)))
         ipar=mod(jtot+j1i+j2i+j12j,2)
         facj=facj*fpar(ipar)
         do il=0,lamax
           ll=l(il)
           ll1=l1(il)
           ll2=l2(il)
           facl=dsqrt(dfloat((2*ll1+1)*(2*ll2+1)))*(2*ll+1)
           som1=threej0(ll,lllj,llli)
           if(dabs(som1)<eps) cycle
           som2=threej0(ll1,j1j,j1i)
           if(dabs(som2)<eps) cycle
           som3=threej0(ll2,j2j,j2i)
           if(dabs(som3)<eps) cycle
           som4=sixj(lllj,llli,ll,j12i,j12j,jtot)  ! I keep this ...
           if(dabs(som4)<eps) cycle
           som=som1*som2*som3*som4*ninej(j12j,j2j,j1j,j12i,j2i,j1i,ll,ll2,ll1)
           fl(i,j,il)=som*faci*facj*facl
         enddo
      enddo
   enddo
   !$OMP END PARALLEL DO


   if(printw) then
     irest=(nlev/10)
     do lam=0,lamax
        do j=1,irest
           write(47,'(2I10)') lam,(j-1)*10+1
           do i=1,nlev
              write(47,'(I3,10E14.6)') i,fl(i,(j-1)*10+1:j*10,lam)
           enddo
        enddo
        if (irest*10 .NE. nlev) then
           write(47,'(2I10)') lam,(j-1)*10+1
                
           do i=1,nlev
              write(47,'(I3,10E14.6)') i,fl(i,irest*10+1:nlev,lam)
           enddo
        endif
     enddo
  endif

end subroutine V03
