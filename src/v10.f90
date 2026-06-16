subroutine V10(jtot,nlev,fl)
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
  implicit real(kind=8) (a-h,o-z)
  integer :: i, j, k,  ndim, lam, irest
  integer :: lmin,lmax
  integer, intent(in) :: jtot,nlev
  real(kind=dp) :: fac,som, f3j0
  real(kind=dp), intent(out), dimension(nlev,nlev,0:lamax) :: fl
  logical :: printw
  real(kind=qp) :: w3js, w6js, w9js, som1,ninej
 ! write(9,*) fact
    printw=.false.
! printw=.true.
  eps=1.d-10

  fl(:,:,:)=0.d0
  do i=1,nlev
        j1i=j1tab(i)
        kai=katab(i)
        j2i=j2tab(i)
        j12i=j12tab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((2*j1i+1)*(2*j2i+1)*(2*llli+1)*(2*j12i+1)))
        ieps=1
!        if(kai.NE.0)ieps=(-1)**j1i
        if(kai<0)ieps=-ieps
              if(kai==0)then
                      ikai=2
              else
                      ikai=1
              endif
              kai=abs(kai)
     do j=i,nlev
        j1j=j1tab(j)
        kaj=katab(j)
        j2j=j2tab(j)
        j12j=j12tab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((2*j1j+1)*(2*j2j+1)*(2*lllj+1)*(2*j12j+1)))
        jeps=1
!        if(kaj.NE.0)jeps=(-1)**j1j
        if(kaj<0)jeps=-jeps
              if(kaj==0)then
                      ikaj=2
              else
                      ikaj=1
              endif
              kaj=abs(kaj)
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
        ll1=l1(il)
        ll2=l2(il)
!        fack=dsqrt(dfloat(2*ll+1))
        fack=dsqrt(dfloat(2*ll2+1))*(dfloat(2*ll+1))!/((4.d0*pi))
!        fack=(dfloat(2*ll+1))/((4.d0*pi))
!        fack=dsqrt(dfloat(2*ll2+1))*(dfloat(2*ll+1)/pi)
!        fack=1.d0
        ipar=mod(jtot+j1j+j2j+j12i+kaj,2)
        fack=fack*fpar(ipar)
        ipar=mod(mm1+j1i+j1j+ll2+ll,2)
        fack=fack*(1.d0+dfloat(ieps*jeps)*fpar(ipar))/2.d0
        if((kaj-kai).GE.0)then
           omega=1.d0
           mm1p=+mm1
        else
           ipar=mod(ll+ll1+ll2+mm1,2)
!           ipar=mod(mm1,2)
           omega=fpar(ipar)
           mm1p=-mm1
        endif

              som1=threej0(ll,llli,lllj)*threej0(ll2,j2i,j2j)
              if(dabs(som1)<eps) cycle
              som2=omega*threej(ll1,j1i,j1j,mm1p,(kai),-(kaj))+dfloat(ieps)*threej(ll1,j1i,j1j,mm1,-(kai),-(kaj)) 
              if(dabs(som2)<eps) cycle
              som3=sixj(lllj,llli,ll,j12i,j12j,jtot)
              if(dabs(som3)<eps) cycle
              som4=ninej(j12i,j2i,j1i,j12j,j2j,j1j,ll,ll2,ll1)
              som=som1*som2*som3*som4

              if(mm1==0)then
                      imm1=1
              else
                      imm1=1
              endif
              facm=dsqrt(dfloat(ikai*ikaj))
              fl(i,j,il)=som*faci*facj*fack/(dfloat(imm1)*dsqrt(facm))
!             fl(j,i,il)=fl(i,j,il)
!              write(9,*)i,j,il,fl(i,j,il),som,faci,facj,facl,facm,omega,ieps
           enddo

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
     write(9,'(5I4)')lam,l(lam),m1(lam),l1(lam),l2(lam)
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
end subroutine V10
