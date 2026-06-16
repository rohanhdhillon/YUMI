subroutine c0140
 use mod_constantes
  use mod_pot
  use mod_base
  use mod_input
  !$ use omp_lib
  integer :: i, j, kk,isym1,isym2
  integer :: nlev,lmin,lmax,jmin,jmax,jtot,nopen,nopen1
  real(kind=dp),allocatable, dimension(:,:,:) :: fl
  real(kind=dp),allocatable, dimension(:,:) :: vvv
  real(kind=dp),allocatable, dimension(:,:) :: RP,Tr,Ti,sigma,sigmao
  integer , allocatable, dimension(:) :: j1tabs, katabs, j2tabs, ka2tabs
  real(kind=dp),allocatable, dimension(:) :: RP0
  real(kind=dp) :: Be, alphae, De, mu,A,B,C,A2,B2,C2
  real(kind=dp) :: deltar, rmin1
  integer :: parite, irest,iprint, imp, output, input8, ipot
  logical :: printw,fla
  real time_begin,time_end
  namelist /basis/jmin,jmax,jstep,Be1,Be2,alphae1,alphae2,&
       De1,De2,mu,iprint,j1min,j1max,j1step,j2min,j2max,k1max,j2step,Bc1,&
       isym1,isym2,A,B,C,A2,B2,C2

!**********************************************************************************
!               Toupie symétrique -  Rotateur Linéaire   
!
!      Rist     et al. J. Chem. Phys. 98, 4662 (1993); doi: 10.1063/1.464970 
!
!**********************************************************************************
!
!**********************************************************************************
!
!   lecture de la base (namelist  base)
!
!**********************************************************************************

  input8=8
  output=9
  imp=output
  ipot=10
  printw=.false.
  jtot=3
  j1min=0
  j1max=3
  j1step=1
  j2min=0
  j2max=3
  j2step=1
  k1max=j1max
  isym1=1
  isym2=1
  Be1=0.d0
  alphae1=0.d0
  De1=0.d0
  Be2=0.d0
  alphae2=0.d0
  De2=0.d0
  A=1.d0
  B=1.d0
  C=1.d0
  A2=1.d0
  B2=1.d0
  C2=1.d0
  iprint=0
  rewind(8)
  read(8,nml=basis)
  Be1=Be1-alphae1/2.d0
  Bc1=Bc1-Be1
  Be2=Be2-alphae2/2.d0
  mu=mu*FMPRT
!**********************************************************************************
!
!   construction de la base
!   détermination de la dimension de la base
!
!**********************************************************************************
     allocate(valp(0:j1max,2*j1max+1),vec(0:j1max,2*j1max+1,2*j1max+1))
     allocate(alpha(0:j1max,2*j1max+1),epsk(0:j1max,2*j1max+1))
     allocate(valp2(0:j2max,2*j2max+1),vec2(0:j2max,2*j2max+1,2*j2max+1))
     allocate(alpha2(0:j2max,2*j2max+1),epsk2(0:j2max,2*j2max+1))
     valp=0.d0
     vec=0.d0
     do j1=j1min,j1max,j1step
        jp1=j1+j1+1
        call asymm(j1,valp(j1,1:jp1),vec(j1,1:jp1,1:jp1),alpha(j1,1:jp1),epsk(j1,1:jp1),A,B,C)
     enddo
     valp2=0.d0
     vec2=0.d0
     do j2=j2min,j2max,j2step
        jp2=j2+j2+1
        call asymm(j2,valp2(j2,1:jp2),vec2(j2,1:jp2,1:jp2),alpha2(j2,1:jp2),epsk2(j2,1:jp2),A2,B2,C2)
     enddo
  jmax=0
  do j2=j2min,j2max,j2step
     jp2=j2+j2+1
     do ik2=1,jp2
     ipar=mod((ik2+j2+1),2)
     if(ipar/=isym2)cycle
     do j1=j1min,j1max,j1step
        jp1=j1+j1+1
     do ik=1,jp1
     ipar=mod((ik+j1+1),2)
     if(ipar==isym1)jmax=jmax+1
!     write(output,*)epsk(j1,ik)
!     if(mod(ik,2)==kstep)jmax=jmax+1
!     if(alpha(j1,ik)*epsk(j1,ik)==dfloat(kstep))jmax=jmax+1
!     if(epsk(j1,ik)==dfloat(kstep))jmax=jmax+1
!             jmax=jmax+1
     enddo
     enddo
     enddo
  enddo
  !**********************************************************************************
  write(output,*) ' j1  ka   kc  tau par  j2  ka  kc  tau par     Ej(cm-1)',jmax,' voies totales'
  allocate(k2tab(jmax))
  allocate(ktab(jmax))
  allocate(jtab(jmax))
  allocate(j1tabs(jmax))
  allocate(j2tabs(jmax))
  allocate(katabs(jmax))
  allocate(ka2tabs(jmax))
  do i=1,jmax
     jtab(i)=i
  enddo
!
  !**********************************************************************************
  !
  !       K**2
  !
  !**********************************************************************************
  !
  bbb=0.d0
  nlevo=0
  jmax=0
  do j2=j2min,j2max,j2step
     jp2=j2+j2+1
     do ik2=1,jp2
     ipar=mod((ik2+j2+1),2)
     if(ipar/=isym2)cycle
     do j1=j1min,j1max,j1step
        jp1=j1+j1+1
     do ik=1,jp1
!     ipar=mod(abs(ik),2)
!     if(ipar==kstep)then
!     if(mod(ik,2)==kstep)then
!     if(alpha(j1,ik)*epsk(j1,ik)==dfloat(kstep))then
!     if(epsk(j1,ik)==dfloat(kstep))then
     ipar=mod((ik+j1+1),2)
     if(ipar==isym1)then
             jmax=jmax+1
        aaa= valp(j1,ik)+valp2(j2,ik2)
        if(jmax==1)bbb=aaa
        write(output,'(10I4,2F18.7)') j1,(ik)/2,(jp1-ik+1)/2,ik,int(epsk(j1,ik)),j2,(ik2)/2,(jp2-ik2+1)/2,ik2,int(epsk(j2,ik2)), aaa
     flush(output)
        k2tab(jmax)=2.d0*mu*(E-aaa+bbb)/uacmm1
     if(k2tab(jmax)>0.d0) nlevo=nlevo+1
     j1tabs(jmax)=j1
     j2tabs(jmax)=j2
     katabs(jmax)=ik
     ka2tabs(jmax)=ik2
     endif
     enddo
     enddo
     enddo
  enddo
  !**********************************************************************************
  !
  !       ktab(i) K**2 trié           jtab(i) numéro du niveau i
  !
  !**********************************************************************************
  !
  write(output,*)
  write(output,*) ' j1  ka  kc  j2  ka  kc         k2(cm-1)',nlevo,' voies ouvertes'
  do i=1,jmax-1
     m=i
     do j=i+1,jmax
        if(k2tab(j).GT.(k2tab(m)))m=j
     enddo
     c=k2tab(i)
     k2tab(i)=k2tab(m)
     k2tab(m)=c
     jj=jtab(i)
     jtab(i)=jtab(m)
     jtab(m)=jj
     jj=j1tabs(i)
     j1tabs(i)=j1tabs(m)
     j1tabs(m)=jj
     jj=j2tabs(i)
     j2tabs(i)=j2tabs(m)
     j2tabs(m)=jj
     jj=katabs(i)
     katabs(i)=katabs(m)
     katabs(m)=jj
     jj=ka2tabs(i)
     ka2tabs(i)=ka2tabs(m)
     ka2tabs(m)=jj
     write(output,'(6I4,F18.7)') j1tabs(i),katabs(i)/2,(2*j1tabs(i)-katabs(i)+2)/2, &
             j2tabs(i),ka2tabs(i)/2,(2*j2tabs(i)-ka2tabs(i)+2)/2,k2tab(i)*uacmm1/(2.d0*mu) 
!     write(output,*)katabs
     flush(output)

  enddo
  i=jmax
     write(output,'(6I4,F18.7)') j1tabs(i),katabs(i)/2,(2*j1tabs(i)-katabs(i)+2)/2, &
             j2tabs(i),ka2tabs(i)/2,(2*j2tabs(i)-ka2tabs(i)+2)/2,k2tab(i)*uacmm1/(2.d0*mu) 
  ktab=k2tab
  deallocate(k2tab)

!**********************************************************************************
!
!   pas d'intégration
!
!**********************************************************************************
  drr=(pi/dsqrt(E*(2.d0*mu)/uacmm1))/dfloat(npas)
  npas=int((rmax-rmin)/drr)
  if((npas/2)*2.NE.npas) then
     npas=npas+1
     drr=((rmax-rmin)/dfloat(npas))
  endif

!**********************************************************************************
!
!   lecture des Vlambda
!
!**********************************************************************************
  read(ipot,*)
  read(ipot,*) lamax,nmax,rmin1,deltar
!  write(output,*) lamax,nmax,rmin1,deltar
!     flush(output)
  allocate(ri(0:nmax),hh(0:nmax))
  allocate(fi(0:nmax,0:lamax),fs(0:nmax,0:lamax))
  allocate(l1(0:lamax),m1(0:lamax),l2(0:lamax),m2(0:lamax),l(0:lamax))
  allocate(vvl(0:lamax))
  allocate(vvv(60,0:lamax))

    irf=2
  select case(irf)
    case(1)

  ri(0)=rmin1!/0.529177d0
  deltar=deltar!/0.529177d0
  do i=1,nmax
     ri(i)=ri(0)+dfloat(i)*deltar
  enddo
    case(2)
  read(ipot,*) ri
!  ri=ri/0.529177d0
  end select

  !  rmax=ri(nmax)
  do i=0,nmax-1
     hh(i)=ri(i+1)-ri(i)
  enddo
  fla=.false.
  if(fla)then
     it=0
     do i=0,l1max,2
        do j=0,l2max
           mmin=iabs(j-i)
           do il=mmin,i+j,2
              l1(it)=i
              l2(it)=j
              l(it)=il
              !        write(*,*)it,i,j,il,l1max,l2max
              it=it+1
           enddo
        enddo
     enddo
!     write(output,*)'it, lamax=',it,lamax
     do il=0,lamax
        READ(ipot,*) fi(0:nmax,il)
        !        write(*,*)il
        !        write(*,*) fi(0:nmax,il)
     enddo

  else
!  write(output,*) lamax,lamax,nmax,rmin1,deltar
     do il=0,lamax
        read(ipot,*) l1(il),m1(il),l2(il),m2(il),l(il)
!        write(output,*) l1(il),m1(il),l2(il),l(il)
           READ(ipot,*) fi(0:nmax,il)
!           write(OUTPUT,*) fi(irest*100:nmax,il)
     enddo
  endif
!**********************************************************************************
!
!   conversion des Vlambda   (cm-1 ------> ua)
!
!**********************************************************************************
  fi=fi*(2.d0*mu)!/uacmm1
!**********************************************************************************
!
!   calcul des dérivés seconde des Vlambda
!
!**********************************************************************************
  do il=0,lamax
     call spline1(ri,fi(:,il),fs(:,il),nmax)
  enddo
!        r=3.0d0
!        do ir=1,60
!        call pot
!        vvv(ir,:)=vvl
!        r=r+0.2d0
!        enddo
!  do il=0,lamax
!        write(output,'(4I5)') l1(il),m1(il),l2(il),l(il)
!        write(output,'(100F19.8)') vvv(:,il)
!  enddo

  allocate(sigma(jmax,jmax))
  allocate(sigmao(nlevo,nlevo))
  sigma=0.d0
  !**********************************************************************************
  !
  !............................................. Boucle sur J total
  !
  !**********************************************************************************
  write(output,*)'                      voies          voies ' 
  write(output,*)'      J  par         totales        ouvertes       Temps CPU'
  flush(output)
  do jtot=jtotl,jtotu,jtstep
     !**********************************************************************************
     !
     !............................................. Boucle sur la parité
     !
     !**********************************************************************************
     j12step=1
     do parite=0,0
        time_begin = OMP_GET_WTIME()
!**********************************************************************************
!
!  détermination du nombre de voies couplées
!
!**********************************************************************************
        nlev=0
        do i=1,jmax   
           j1=j1tabs(i)
           j2=j2tabs(i)
           ka=katabs(i)
           ka2=ka2tabs(i)
           j12min=iabs(j1-j2)
           j12max=j1+j2
           do j12=j12min,j12max,j12step
              lmin=iabs(jtot-j12)
              lmax=jtot+j12
              do il=lmin,lmax
!                 kk=jtot+j1+ka+il+parite
!                 kk=j1+ka+parite
!kk=0
!ika=mod(abs(ka),2)
                 kk=jtot+j1+j2+il+parite
!                 if(int(alpha(j1,ka)*epsk(j1,ka))>0)kk=kk+1
!                 if(epsk(j1,ka)>0)kk=kk+1
!                 if(alpha(j1,ka)>0)kk=kk+1
                 kk=iabs(kk)
                 kk=0
                    if(mod(kk,2).eq.0) then
                    nlev=nlev+1
                    endif
              end do
           enddo
        enddo
        if (nlev/=0) then
           allocate(j1tab(nlev))
           allocate(katab(nlev))
           allocate(ka2tab(nlev))
           allocate(j2tab(nlev))
           allocate(j12tab(nlev))
           allocate(jjtab(nlev))
           allocate(ltab(nlev))
           allocate(k2tab(nlev))
           allocate(fl(nlev,nlev,0:lamax))
           allocate(RP(nlev,nlev))
           allocate(RP0(nlev))

!**********************************************************************************
!
!   détermination des voies couplées
!
!**********************************************************************************

           nlev=0
           do i=1,jmax   
           j1=j1tabs(i)
           j2=j2tabs(i)
           ka=katabs(i)
           ka2=ka2tabs(i)
              j12min=iabs(j1-j2)
              j12max=j1+j2
              do j12=j12min,j12max,j12step
                 lmin=iabs(jtot-j12)
                 lmax=jtot+j12
                 do il=lmin,lmax
                 kk=jtot+j1+j2+il+parite
!                 if(int(alpha(j1,ka)*epsk(j1,ka))>0)kk=kk+1
!                 if(epsk(j1,ka)>0)kk=kk+1
!                 if(alpha(j1,ka)>0)kk=kk+1
                 kk=iabs(kk)
                 kk=0
                    if(mod(kk,2).eq.0) then
                       nlev=nlev+1
                       j1tab(nlev)=j1
                       katab(nlev)=ka
                       ka2tab(nlev)=ka2
                       j2tab(nlev)=j2
                       j12tab(nlev)=j12
                       ltab(nlev)=il
                       k2tab(nlev)=ktab(i)
                       jjtab(nlev)=jtab(i)
                    endif
                 enddo
              enddo
           enddo
!**********************************************************************************
!
!   calcul de la matrice de couplage  fl
!
!**********************************************************************************
           fl(:,:,:)=0.d0
           call V140(jtot,nlev,fl)


!**********************************************************************************
!
!   propagateur
!
!**********************************************************************************
                call ijohnson(rmin,rmax,RP,nlev,npas,fl)
           do i=1,nlev
              do j=i,nlev
                 RP(i,j)=(RP(i,j)+RP(j,i))/2.d0
                 RP(j,i)=RP(i,j)
              enddo
           enddo
!**********************************************************************************
!
!  détermination du nombre de voies couplées ouvertes
!
!**********************************************************************************
           nopen=0
           do j=1,nlev
              if(k2tab(j)>0.d0) then
                 nopen=nopen+1
              endif
           enddo




           if(nopen.EQ.0) then
              deallocate(ltab,fl,RP,k2tab,j1tab,katab,ka2tab,j2tab,j12tab,jjtab,RP0)
              cycle
           endif


!**********************************************************************************
!
!   calcul de la matrice    T  = Tr  +  i Ti
!
!**********************************************************************************
           nopen1=nopen
           allocate(Tr(nopen1,nopen1),Ti(nopen1,nopen1))
           call tmat(RP,nlev,nopen1,ltab,k2tab,rmax,Tr,Ti,iprint,imp)
!**********************************************************************************
!
!   calcul des sections efficaces
!
!**********************************************************************************
           Tr=Tr*Tr
           Ti=Ti*Ti
           Tr=(Tr+Ti)*dfloat(2*jtot+1)*dfloat(jtstep)*pi
           do i=1,nopen
              do j=1,nopen
           j1=j1tab(i)
           j2=j2tab(i)
           ka=katab(i)
                 sigma(jjtab(i),jjtab(j))=sigma(jjtab(i),jjtab(j))+Tr(i,j)/(k2tab(i)*dfloat((2*j1+1)*(2*j2+1)))
              enddo
           enddo
           deallocate(ltab,fl,RP,Tr,Ti,k2tab,j1tab,katab,ka2tab,j2tab,j12tab,jjtab,RP0)
        else
           nopen=0
        endif
  time_end = OMP_GET_WTIME()
  write(output,'(I8,I3,2I15,F19.3)') jtot,parite,nlev,nopen,time_end-time_begin
  flush(output)
     enddo
     !**********************************************************************************
     !
     !..........................................     fin parité
     !
     !**********************************************************************************
  enddo
  !**********************************************************************************
  !
  !..........................................     fin J total
  !
  !**********************************************************************************
!**********************************************************************************
!
!   conversion Angsetrum**2 <--- ua**2
!
!**********************************************************************************
  sigma=sigma*A02A2
  write(output,*)
!  write(output,*)'      N°    (  J1,  ka,  J2)          EJ(cm-1)'
  write(output,*) ' N°  j2   j1  ka  kc        k2(cm-1)'
do i=1,nlevo
           m=jtab(i)
           j1=j1tabs(i)
           j2=j2tabs(i)
           ka=katabs(i)
     write(output,'(5I4,F18.7)') i,j2,j1,ka/2,(2*j1-ka+2)/2,ktab(i)*uacmm1/(2.d0*mu) 
!write(output,2500) i,m,j1,ka,j2, ktab(i)*uacmm1/(2.d0*mu)
enddo
  sigmao(:,:)=sigma(jtab(1:nlevo),jtab(1:nlevo))
  write(output,*)
  call imprime2(iprint,imp,nlevo,sigmao,'sigma     ')
2500   format(2I5,'   (',I4,',',I4,',',I4,')   ',F16.4)
end subroutine c0140
