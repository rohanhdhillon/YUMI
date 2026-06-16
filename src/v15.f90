subroutine V15(jtot,nlev,fl)
!$  use OMP_LIB
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
  implicit real(kind=8) (a-h,o-z)
  integer :: i, j, k, lam
  integer, intent(in) :: jtot,nlev
  real(kind=dp) :: fac,som, f3j0
  real(kind=dp), intent(out), dimension(nlev,nlev,0:lamax) :: fl
  real(kind=qp) :: w3js, w6js, w9js, som1,som2,som3,som4,som5
  logical :: printw
 printw=.false.
  eps=1.d-10
 sqr2=dsqrt(0.5d0)
  fl(:,:,:)=0.d0
  do i=1,nlev
        j1i=j1tab(i)
        kai=katab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((2*j1i+1)*(2*llli+1)))
!        faci=dsqrt(dfloat((2*j1i+1)*(2*llli+1)/(pi)))
        ieps=ptab(i)
!        write(9,*)j,ieps
        aieps=fpar(ieps)
!        aieps=1.d0
!        aieps=sign(aieps,kai)!*(-1)**j1j
        if(kai<0)aieps=-aieps
        kai=abs(katab(i))
!        write(9,*)j,aieps
              if(kai==0)then
                      facmi=sqr2
              else
                      facmi=1.d0
              endif
     do j=i,nlev
        j1j=j1tab(j)
        kaj=katab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((2*j1j+1)*(2*lllj+1)))
        jeps=ptab(j)
        ajeps=fpar(jeps)
!        ajeps=1.d0
!        ajeps=sign(ajeps,kaj)!*(-1)**j1j
        if(kaj<0)ajeps=-ajeps
        kaj=abs(katab(j))
              if(kaj==0)then
                      facmj=sqr2
              else
                      facmj=1.d0
              endif
!
!
!*!$OMP PARALLEL DO PRIVATE(il, ll, mm1, ipar, fack, som1, som6, som5,omega,mm1p,facmm)
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
              if(mm1==0)then
                      facmm=sqr2
              else
                      facmm=1.d0
              endif
        ipar=mod(iabs(-jtot+j1i+j1j+kaj),2)
        fack=fpar(ipar)
        ipar=mod(mm1+j1i+j1j+ll,2)
        fack=fack*(1.d0+aieps*ajeps*fpar(ipar))/2.d0
        if((kaj-kai).GE.0)then
           omega=1.d0
           mm1p=-mm1
        else
           ipar=mod(mm1,2)
           omega=fpar(ipar)
           mm1p=+mm1
        endif

              som5=threej0(llli,lllj,ll)
              if(dabs(som5)<eps) cycle
              som1=omega*threej(j1i,j1j,ll,kai,-kaj,mm1p)+aieps*threej(j1i,j1j,ll,-kai,-kaj,mm1)
              if(dabs(som1)<eps) cycle
              som6=(sixj(j1j,lllj,jtot,llli,j1i,ll))
              if(dabs(som6)<eps) cycle
              fl(i,j,il)=som1*som5*som6*faci*facj*fack*dsqrt(dfloat(2*ll+1)/(4.d0*pi))*facmi*facmj*facmm
           enddo
!*!$OMP END PARALLEL DO
!
!
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

end subroutine V15
