subroutine V4(jtot,nlev,fl)
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
  real(kind=qp) :: w3js, w6js, w9js, som1
 ! write(9,*) fact
 !   printw=.true. 
 printw=.true.
 aa1=pi**(1.5d0)
 aa=pi*dsqrt(pi)
 !pi15=pi**(1.5d0)
 aa2=5.56832799683171d0

  fl(:,:,:)=0.d0
  do i=1,nlev
        j1i=j1tab(i)*2
        kai=katab(i)*2
        j2i=j2tab(i)*2
        j12i=j12tab(i)*2
        llli=ltab(i)*2
        faci=dfloat((j1i+1)*(j2i+1))*dsqrt(dfloat((llli+1)))
        ieps=1
        ieps=isign(ieps,kai)*(-1)**j1i
!        write(9,*)j1i,kai,j2i,j12i,llli
              if(kai==0)then
                      ikai=2
              else
                      ikai=1
              endif
     do j=i,nlev
        j1j=j1tab(j)*2
        kaj=katab(j)*2
        j2j=j2tab(j)*2
        j12j=j12tab(j)*2
        lllj=ltab(j)*2
        facj=dfloat((j1j+1)*(j2j+1))*dsqrt(dfloat((lllj+1)))
        jeps=1
        jeps=isign(jeps,kaj)*(-1)**j1j
!        write(9,*)j1i,kai,j2i,j12i,llli
!        write(9,*)ieps,jeps
              if(kaj==0)then
                      ikaj=2
              else
                      ikaj=1
              endif
        do il=0,lamax
        ll=l(il)*2
        mm1=m1(il)*2
        ll1=l1(il)*2
        ll2=l2(il)*2
        fack=(-1.d0)**(jtot+(kaj+mm1)/2)
!        fack=fack*(1.d0+dfloat(ieps*jeps)*((-1.d0)**((mm1+j1i+j1j+ll2+ll)/2)))
              som1=w3js(llli,ll,lllj,0,0,0)*w3js(j2i,ll2,j2j,0,0,0)
              som1=som1*(w3js(j1i,ll1,j1j,kai,-mm1,-kaj)+((-1.d0)**(mm1/2))*w3js(j1i,ll1,j1j,kai,mm1,-kaj)) 
              som1=(som1*w6js(j12j,lllj,2*jtot,llli,j12i,ll))
              som1=(som1*w9js(ll1,j1j,j1i,ll2,j2j,j2i,ll,j12j,j12i))
              som=real(som1,8)
!              fl(i,j,il)=som*faci*facj*facl/(pi**1.5d0) 
              if(mm1==0)then
                      imm1=2
              else
                      imm1=1
              endif
              facm=dfloat(imm1)
              fl(i,j,il)=som*faci*facj*fack/((facm)) 
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
end subroutine V4
