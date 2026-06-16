subroutine V04sp(jtot,nlev,nsp)
!$ use OMP_LIB
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_fact
  implicit real(kind=8) (a-h,o-z)
  integer :: i, j, k,  ndim, lam, irest, taui,tauj,kai,kaj
  integer :: j1i,j1j,j2i,j2j,j12i,j12j,llli,lllj,ll,mm1,ll1,ll2
  integer, intent(in) :: jtot,nlev
  integer, intent(out) :: nsp
  real(kind=dp) :: fac,som,som1,som2,som3,som4, f3j0, somt
  real(kind=dp), dimension(nlev,nlev,0:lamax) :: fl
  integer, allocatable, dimension(:) :: indexl2,indexj2
!  real(kind=dp), allocatable, intent(out) :: flsp(:)
!  integer, allocatable, intent(out) :: indexl(:), indexj(:)
  logical :: printw
  real(kind=qp) :: w3js, w6js, w9js, ninej,sixj,threej0,threej
 ! write(9,*) fact
    printw=.false.
!    printw=.true.
  eps=1.d-10
  fpim1=1.d0/(4.d0*pi)

!  nbv2=(nlev*nlev*(lamax+1))/6
!  allocate(fl(nbv2),indexl2(nbv2),indexj2(nbv2))
  do i=1,nlev
        j1i=j1tab(i)
        taui=katab(i)
        j2i=j2tab(i)
        j12i=j12tab(i)
        llli=ltab(i)
        faci=dsqrt(dfloat((j1i+j1i+1)*(j2i+j2i+1)*(llli+llli+1)*(j12i+j12i+1)))
                kaimax=j1i
        if((alpha(j1i,taui)<0.and.mod(j1i,2)/=0).or.alpha(j1i,taui)>0.and.mod(j1i,2)==0)kaimax=j1i-1
!     do j=i,nlev
     do j=1,i
        j1j=j1tab(j)
        tauj=katab(j)
        j2j=j2tab(j)
        j12j=j12tab(j)
        lllj=ltab(j)
        facj=dsqrt(dfloat((j1j+j1j+1)*(j2j+j2j+1)*(lllj+lllj+1)*(j12j+j12j+1)))
                kajmax=j1j
        if((alpha(j1j,tauj)<0.and.mod(j1j,2)/=0).or.alpha(j1j,tauj)>0.and.mod(j1j,2)==0)kajmax=j1j-1
!$OMP PARALLEL
! rien!$OMP FIRSTPRIVATE(l,l1,l2,m1)
!$OMP DO SCHEDULE(DYNAMIC,2) PRIVATE(il,ll,mm1,ll1,ll2,fack,som1,som3,som4,som2,somt,som,kai,kaj,ipar1,ipar)
        do il=0,lamax
        ll=l(il)
        mm1=m1(il)
        ll1=l1(il)
        ll2=l2(il)
        if(mod(llli+lllj+ll,2)/=0)cycle
        if(mod(j2i+j2j+ll2,2)/=0)cycle
        if(ll<abs(llli-lllj).or.ll>abs(llli+lllj)) cycle
        if(ll2<abs(j2i-j2j).or.ll2>abs(j2i+j2j)) cycle
        if(ll1<abs(j1i-j1j).or.ll1>abs(j1i+j1j)) cycle
!        if(ll<abs(j12i-j12j).or.ll>abs(j12i+j12j)) cycle
        fack=dsqrt(dfloat(ll2+ll2+1)*(dfloat(ll+ll+1)))
              som1=threej0(llli,lllj,ll)*threej0(j2i,ll2,j2j)
!              if(dabs(som1)<eps) cycle
              som3=sixj(llli,lllj,ll,j12j,j12i,jtot)
              if(dabs(som3)<eps) cycle
              som4=ninej(j12j,j12i,ll,j1j,j1i,ll1,j2j,j2i,ll2)
              if(dabs(som4)<eps) cycle
        ipar1=mod(mm1+ll1+ll2+ll,2)
              som=0.d0
        do kai=-kaimax,kaimax,2
        do kaj=-kajmax,kajmax,2
        if(kai/=mm1+kaj.and.kai/=-mm1+kaj)cycle
        ipar=mod(abs(jtot-j1j+j2j-j12i+kai-ll),2)

              som2=threej(j1i,ll1,j1j,-(kai),mm1,(kaj))+fpar(ipar1)*threej(j1i,ll1,j1j,-(kai),-mm1,(kaj)) 
!              if(dabs(som2)<eps) cycle
              som=som+som2*fpar(ipar)*(vec(j1i,kai+j1i+1,taui)*vec(j1j,kaj+j1j+1,tauj))
         enddo
         enddo
              if(mm1==0)som=som*0.5d0
              somt=som*som1*som3*som4*faci*facj*fack*fpim1
              if(dabs(somt)<eps)cycle
              fl(i,j,il)=somt
         enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

     enddo
  enddo
nsp=0
do i=1,nlev
do j=1,i
do il=0,lamax
if(dabs(fl(i,j,il))>eps)nsp=nsp+1
enddo
enddo
enddo
allocate(flsp(nsp),indexl(nsp),indexj(nsp))
ij=-1
ii=0
do i=1,nlev
do j=1,i
ij=ij+1
do il=0,lamax
if(dabs(fl(i,j,il))>eps)then
        ii=ii+1
        flsp(ii)=fl(i,j,il)
        indexl(ii)=il
        indexj(ii)=ij
endif
enddo
enddo
enddo
!write(9,*)nsp
!write(9,*)flsp
!write(9,*)indexl
!write(9,*)indexj
end subroutine V04sp
