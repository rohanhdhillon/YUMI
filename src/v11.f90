subroutine V11(jtot,nlev,fl)
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
  real(kind=dp),intent(out),dimension(nlev,nlev,0:lamax)::fl
  logical :: printw
  real(kind=qp) :: w3js, w6js, w9js, som1,ninej
 ! write(9,*) fact
 !   printw=.true. 
 printw=.false.
 aa1=pi**(1.5d0)
 aa=pi*dsqrt(pi)
 !pi15=pi**(1.5d0)
 aa2=5.56832799683171d0

  fl(:,:,:)=0.d0
  do i=1,nlev
        j1i=j1tab(i)
        kai=katab(i)
        j2i=j2tab(i)
        j12i=j12tab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((2*j1i+1)*(2*j2i+1)*(2*llli+1)*(2*j12i+1)))
        ieps=1
        ieps=isign(ieps,kai)!*(-1)**(j1i/2)
        kai=abs(kai)
!        write(9,*)j1i,kai,j2i,j12i,llli
              if(kai==0)then
                      ikai=2
              else
                      ikai=1
              endif
     do j=i,nlev
        j1j=j1tab(j)
        kaj=katab(j)
        j2j=j2tab(j)
        j12j=j12tab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((2*j1j+1)*(2*j2j+1)*(2*lllj+1)*(2*j12j+1)))
        jeps=1
        jeps=isign(jeps,kaj)!*(-1)**(j1j/2)
        kaj=abs(kaj)
!        write(9,*)j1i,kai,j2i,j12i,llli
!        write(9,*)ieps,jeps
              if(kaj==0)then
                      ikaj=2
              else
                      ikaj=1
              endif
!$OMP PARALLEL DO PRIVATE (il,ll,mm1,ll1,ll2,fack,som1,som2,som3,som4,som,facm,imm1)
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
        ll1=l1(il)
        ll2=l2(il)
        fack=dsqrt(dfloat(2*ll+1))
        ipar=mod(jtot+ll1-ll2+j1j-j2j+j12i-lllj-llli-kaj-mm1,2)
        fack=fack*fpar(ipar)
        ipar=mod(mm1+j1i+j1j+ll2+ll,2)
        fack=fack*(1.d0+dfloat(ieps*jeps)*fpar(ipar))
              som1=threej0(llli,ll,lllj)*threej0(j2i,ll2,j2j)
              if(dabs(som1)<eps) cycle
              som2=threej(j1i,ll1,j1j,-(kai),mm1,(kaj))+dfloat(ieps*jeps)*threej(j1i,ll1,j1j,(kai),mm1,-(kaj)) &
              +dfloat(ieps)*threej(j1i,ll1,j1j,(kai),mm1,(kaj))+dfloat(jeps)*threej(j1i,ll1,j1j,-(kai),mm1,-(kaj))
              if(dabs(som2)<eps) cycle
              som3=sixj(j12j,lllj,jtot,llli,j12i,ll)
              if(dabs(som3)<eps) cycle
              som4=ninej(j1j,j2j,j12j,j1i,j2i,j12i,ll1,ll2,ll)
              som=som1*som2*som3*som4

              if(mm1==0)then
                      imm1=2
              else
                      imm1=1
              endif
              facm=dfloat(ikai*ikaj)
              fl(i,j,il)=som*faci*facj*fack/(imm1*dsqrt(facm)) 
!             fl(j,i,il)=fl(i,j,il)
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
           do i=1,nlev
              write(9,'(I3,10E14.6)') i,fl(i,(j-1)*10+1:j*10,lam)
           enddo
        enddo
        if (irest*10 .NE. nlev) then
           do i=1,nlev
              write(9,'(I3,10E14.6)') i,fl(i,irest*10+1:nlev,lam)
           enddo
        endif
     enddo
  endif
end subroutine V11
