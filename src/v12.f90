subroutine V12(jtot,nlev,fl)
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
!$ use OMP_LIB  
  implicit real(kind=8) (a-h,o-z)
  integer :: i, j, k,  ndim, lam, irest
  integer :: lmin,lmax
  integer, intent(in) :: jtot,nlev
  real(kind=dp) :: fac,som, f3j0
  real(kind=dp), intent(out), dimension(nlev,nlev,0:lamax) :: fl
  logical :: printw
  real(kind=qp) :: w3js, w6js, w9js, som1
 ! write(9,*) fact
 !   printw=.true. 
 printw=.false.
 aa1=pi**(1.5d0)
 aa=pi*dsqrt(pi)
 !pi15=pi**(1.5d0)
 aa2=5.56832799683171d0
 sqr2=dsqrt(0.5d0)
 eps=1.d-10

  fl(:,:,:)=0.d0
  do i=1,nlev
        j1i=j1tab(i)
        kai=katab(i)
        j2i=j2tab(i)
        j12i=j12tab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((2*j1i+1)*(2*j2i+1)*(2*llli+1)*(2*j12i+1)))
        ipar=ptab(i)
        aieps=fpar(ipar)
!        write(9,*)j1i,kai,j2i,j12i,llli
              if(kai==0)then
                      akai=sqr2
              else
                      akai=1.d0
              endif
     do j=i,nlev
        j1j=j1tab(j)
        kaj=katab(j)
        j2j=j2tab(j)
        j12j=j12tab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((2*j1j+1)*(2*j2j+1)*(2*lllj+1)*(2*j12j+1)))
        ipar=ptab(j)
        ajeps=fpar(ipar)
!        write(9,*)j1i,kai,j2i,j12i,llli
!        write(9,*)ieps,jeps
              if(kaj==0)then
                      akaj=sqr2
              else
                      akaj=1.d0
              endif
!$OMP PARALLEL DO PRIVATE(il,ll,mm1,ll1,ll2,fack,ipar,omega,mm1p,som1,som2,som3,som4,som)              
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
        ll1=l1(il)
        ll2=l2(il)
        fack=dsqrt(dfloat(2*ll2+1))
!        fack=1.d0/dsqrt(2.d0)
!        fack=1.d0
        ipar=(-jtot+j1j+j2j-j12i-kaj)
        ipar=mod(ipar,2)
        fack=fack*fpar(ipar)
        fack=fack*(dfloat(2*ll+1)/(4.d0*pi))
        if((kaj-kai).GE.0)then
           omega=1.d0
           mm1p=-mm1
        else
!           ipar=mod((mm1)/2,2)
           ipar=mod(ll1+ll2+ll+mm1,2)
!           ipar=mod(iabs(kaj-kai)/2,2)
           omega=fpar(ipar)
           mm1p=+mm1
        endif
        ipar=mod(mm1+j1i+j1j-2*kaj+ll+ll2,2)
!        ipar=mod((mm1)/2,2)
        fack=fack*(1.d0+aieps*ajeps*fpar(ipar))/2.d0
!              som1=(omega+dfloat(ieps))*w3js(llli,lllj,ll,0,0,0)
!              if(som1<eps)cycle
!              som2=(omega*w3js(j1i,j1j,ll,kai,-kaj,mm1)+aieps*w3js(j1i,j1j,ll,kai,-kaj,mm1))
!              if(som2<eps)cycle
!              som3=(omega+dfloat(ieps))*(w6js(j1j,lllj,2*jtot,llli,j1i,ll))
!              write(9,*)omega,kai,kaj,mm1/2
!              som=real(som1,8)
!              som=som1*som2*som3
!              fl(i,j,il)=som*faci*facj*facl/(pi**1.5d0) 
!              som1=w3js(lllj,ll,llli,0,0,0)*w3js(j2j,ll2,j2i,0,0,0)
!              if(dabs(som1)<eps) cycle
!              som3=w6js(lllj,llli,ll,j12i,j12j,2*jtot)
!              if(dabs(som3)<eps) cycle
!              som4=w9js(j12i,j2i,j1i,j12j,j2j,j1j,ll,ll2,ll1)
!              if(dabs(som4)<eps) cycle
!              som2=omega*w3js(j1j,ll1,j1i,-kaj,mm1p,kai) &
!                  +aieps*w3js(j1j,ll1,j1i,-kaj,mm1,-kai)
              som1=threej0(lllj,ll,llli)*threej0(j2j,ll2,j2i)
              if(dabs(som1)<eps) cycle
              som2=omega*threej(j1j,ll1,j1i,-kaj,mm1p,kai) &
                  +aieps*threej(j1j,ll1,j1i,-kaj,mm1,-kai)
              if(dabs(som2)<eps) cycle
              som3=sixj(lllj,llli,ll,j12i,j12j,jtot)
              if(dabs(som3)<eps) cycle
              som4=ninej(j12i,j2i,j1i,j12j,j2j,j1j,ll,ll2,ll1)
              som=som1*som2*som3*som4
!              if(mm1==0)then
!                      facmm=0.5d0
!              else
!                      facmm=1.d0
!              endif
             
              fl(i,j,il)=som*faci*facj*fack*akai*akaj
!              write(9,*)i,j,il,fl(i,j,il),som,faci,facj,facl,facm,omega,ieps
           enddo
!$OMP END PARALLEL DO
!              write(9,*)i,j,fl(i,j,0)
     enddo
  enddo
    do lam=0,lamax
       do i=1,nlev
          do j=i+1,nlev
             fl(j,i,lam)=fl(i,j,lam)
          enddo
       enddo
    enddo
  if(printw) then
     irest=(nlev/10)
     do lam=0,lamax
        do j=1,irest
!           print *, (j-1)*10+1
           do i=1,nlev
              write(9,'(I3,10E14.6)') i,fl(i,(j-1)*10+1:j*10,lam)
           enddo
        enddo
        if (irest*10 .NE. nlev) then
!           print *, (j-1)*10+1
           do i=1,nlev
              write(9,'(I3,10E14.6)') i,fl(i,irest*10+1:nlev,lam)
           enddo
        endif
     enddo
  endif
end subroutine V12
