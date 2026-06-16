subroutine V40(jtot,nlev,fl)
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
  implicit real(kind=8) (a-h,o-z)
  integer :: i, j, k,  ndim, lam, irest, taui,tauj,kai,kaj,tau2i,tau2j,ka2i,ka2j
  integer :: j1i,j1j,j2i,j2j,j12i,j12j,llli,lllj,ll,mm1,ll1,ll2,mm2
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
        tau2i=ka2tab(i)
        j12i=j12tab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((j1i+j1i+1)*(j2i+j2i+1)*(llli+llli+1)*(j12i+j12i+1)))
     do j=i,nlev
        j1j=j1tab(j)
        tauj=katab(j)
        j2j=j2tab(j)
        tau2j=ka2tab(j)
        j12j=j12tab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((j1j+j1j+1)*(j2j+j2j+1)*(lllj+lllj+1)*(j12j+j12j+1)))
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
        ll1=l1(il)
        mm2=m2(il)
        ll2=l2(il)
        ipar7=mod(abs(ll1+mm1+ll2+mm2+ll),2)
        if(ll<abs(llli-lllj).or.ll>abs(llli+lllj))cycle
              som1=threej0(llli,lllj,ll)
              if(dabs(som1)<eps) cycle
              som3=sixj(llli,lllj,ll,j12j,j12i,jtot)
              if(dabs(som3)<eps) cycle
        ipar5=mod(abs(ll1+ll2),2)
              som4=ninej(j12j,j12i,ll,j1j,j1i,ll1,j2j,j2i,ll2)
              som5=0.d0
              if(ll1/=ll2)som5=fpar(ipar5)*ninej(j12j,j12i,ll,j1j,j1i,ll2,j2j,j2i,ll1)
!              if(dabs(som4)<eps) cycle
              som=0.d0
        do kai=-j1i,j1i
        if(vec(j1i,kai+j1i+1,taui)==0.d0)cycle
        do kaj=-j1j,j1j
        if(vec(j1j,kaj+j1j+1,tauj)==0.d0)cycle
        do ka2i=-j2i,j2i
        if(vec2(j2i,ka2i+j2i+1,tau2i)==0.d0)cycle
        do ka2j=-j2j,j2j
        if(vec2(j2j,ka2j+j2j+1,tau2j)==0.d0)cycle
        ipar=mod(abs(jtot+j1j+j2j-j12i+kai+ka2i-ll),2)

              som2=threej(j1i,ll1,j1j,-(kai),mm1,(kaj))*threej(j2i,ll2,j2j,-(ka2i),mm2,(ka2j)) 
              som6=0.d0
              som8=0.d0
              if(ll1/=ll2)then
                      som6=threej(j1i,ll2,j1j,-(kai),mm2,(kaj))*threej(j2i,ll1,j2j,-(ka2i),mm1,(ka2j)) 
                      som8=fpar(ipar7)*threej(j1i,ll2,j1j,-(kai),-mm2,(kaj))*threej(j2i,ll1,j2j,-(ka2i),-mm1,(ka2j)) 
              endif
              som7=fpar(ipar7)*threej(j1i,ll1,j1j,-(kai),-mm1,(kaj))*threej(j2i,ll2,j2j,-(ka2i),-mm2,(ka2j)) 
              if(ll1/=ll2)som8=fpar(ipar7)*threej(j1i,ll2,j1j,-(kai),-mm2,(kaj))*threej(j2i,ll1,j2j,-(ka2i),-mm1,(ka2j)) 
!              if(dabs(som2)<eps) cycle
  som=som+((som2+som7)*som4+(som6+som8)*som5)*fpar(ipar)*vec(j1i,kai+j1i+1,taui)*vec(j1j,kaj+j1j+1,tauj) &
          *vec2(j2i,ka2i+j2i+1,tau2i)*vec2(j2j,ka2j+j2j+1,tau2j)

         enddo
         enddo
         enddo
         enddo
              fl(i,j,il)=som*som1*som3*faci*facj
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
end subroutine V40
