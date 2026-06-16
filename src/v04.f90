subroutine V04(jtot,nlev,fl)
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
  implicit real(kind=8) (a-h,o-z)
  integer :: i, j, k,  ndim, lam, irest, taui,tauj,kai,kaj
  integer :: j1i,j1j,j2i,j2j,j12i,j12j,llli,lllj,ll,mm1,ll1,ll2
  integer, intent(in) :: jtot,nlev
  real(kind=dp) :: fac,som,som1,som2,som3,som4, f3j0
  real(kind=dp), intent(out), dimension(nlev,nlev,0:lamax) :: fl
  logical :: printw
  real(kind=qp) :: w3js, w6js, w9js, ninej,sixj,threej0,threej
 ! write(9,*) fact
    printw=.false.
!    printw=.true.
  eps=1.d-10

  fl(:,:,:)=0.d0
  do i=1,nlev
        j1i=j1tab(i)
        taui=katab(i)
        j2i=j2tab(i)
        j12i=j12tab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((j1i+j1i+1)*(j2i+j2i+1)*(llli+llli+1)*(j12i+j12i+1)))
     do j=i,nlev
        j1j=j1tab(j)
        tauj=katab(j)
        j2j=j2tab(j)
        j12j=j12tab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((j1j+j1j+1)*(j2j+j2j+1)*(lllj+lllj+1)*(j12j+j12j+1)))
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
        ll1=l1(il)
        ll2=l2(il)
        if(ll<abs(llli-lllj).or.ll>abs(llli+lllj)) cycle
        if(ll2<abs(j2i-j2j).or.ll2>abs(j2i+j2j)) cycle
        fack=dsqrt(dfloat(ll2+ll2+1)*(dfloat(ll+ll+1)))/((4.d0*pi))
              som1=threej0(llli,lllj,ll)*threej0(j2i,ll2,j2j)
              if(dabs(som1)<eps) cycle
              som3=sixj(llli,lllj,ll,j12j,j12i,jtot)
              if(dabs(som3)<eps) cycle
              som4=ninej(j12j,j12i,ll,j1j,j1i,ll1,j2j,j2i,ll2)
              if(dabs(som4)<eps) cycle
        ipar1=mod(mm1+ll1+ll2+ll,2)
              som=0.d0
        do kai=-j1i,j1i
!        do kai=0,j1i
        if(vec(j1i,kai+j1i+1,taui)==0.d0)cycle
        do kaj=-j1j,j1j
!        do kaj=0,j1j
        if(kai/=mm1+kaj.and.kai/=-mm1+kaj)cycle
        if(vec(j1j,kaj+j1j+1,tauj)==0.d0)cycle
        ipar=mod(abs(jtot-j1j+j2j-j12i+kai-ll),2)

!        write(9,*)
!        write(9,*)'(',j1i,ll1,j1j,')'
!        write(9,*)'(',kai,mm1,kaj,')'
!        write(9,*)ll1,mm1,ll2,ll
!        write(9,*)
!              write(9,*)threej(j1i,ll1,j1j,-(kai),mm1,(kaj)),threej(j1i,ll1,j1j,-(kai),-mm1,(kaj)) 
              som2=threej(j1i,ll1,j1j,-(kai),mm1,(kaj))+fpar(ipar1)*threej(j1i,ll1,j1j,-(kai),-mm1,(kaj)) 
              if(dabs(som2)<eps) cycle
!              som=som+som2*fack*vec(j1i,kai+j1i+1,taui)*vec(j1j,kaj+j1j+1,tauj)
!              write(9,*)som,kai,kaj
              som=som+som2*fpar(ipar)*(vec(j1i,kai+j1i+1,taui)*vec(j1j,kaj+j1j+1,tauj))
!              write(9,*)j1i,kai+j1i+1,taui,kai,som2
!              write(9,*)j1j,kaj+j1j+1,tauj,kaj,som
!              write(9,*)vec(j1i,kai+j1i+1,taui),vec(j1j,kaj+j1j+1,tauj)

         enddo
         enddo
!        write(9,*)'**********************************'
              if(mm1==0)som=som*0.5
!              fl(i,j,il)=fl(i,j,il)+som*som1*som3*som4*faci*facj/(dfloat(imm1))
              fl(i,j,il)=som*som1*som3*som4*faci*facj*fack
         enddo

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
  do i=1,nlev
     write(9,'(6I6)')i,j1tab(i),j2tab(i),j12tab(i),ltab(i),katab(i)
  enddo
     irest=(nlev/10)
     do lam=0,lamax
     write(9,'(5I4)')lam,l1(lam),m1(lam),l2(lam),l(lam)
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
end subroutine V04
