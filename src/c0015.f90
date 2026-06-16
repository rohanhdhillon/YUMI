subroutine c0015
 use mod_constantes
  use mod_pot
  use mod_base
  use mod_input
  !$ use omp_lib
  integer :: i, j, kk
  integer :: nlev,lmin,lmax,jmin,jmax,jtot,nopen,nopen1
  real(kind=dp),allocatable, dimension(:,:,:) :: fl
  real(kind=dp),allocatable, dimension(:,:) :: RP,Tr,Ti,sigma,sigmao
  integer , allocatable, dimension(:) :: j1tabs, katabs
  real(kind=dp) :: Be, alphae, De, mu
  real(kind=dp) :: deltar, rmin1
  integer :: parite, irest,iprint, imp, output, input8, ipot
  logical :: printw,fla
  real time_begin,time_end
  namelist /basis/jmin,jmax,jstep,Be1,Be2,alphae1,alphae2,&
       De1,De2,mu,iprint,j1min,j1max,j1step,j2min,j2max,j2step,Bc1,&
       kmax,kstep

  !**********************************************************************************
  !
  !
  !**********************************************************************************
  !
  input8=8
  output=9
  imp=output
  ipot=10
  printw=.false.
!**********************************************************************************
!
!   lecture de la base (namelist  base)
!
!**********************************************************************************
  
  jtot=3
  j1min=0
  j1max=3
  j1step=1
  kmax=j1max
  kstep=1
  Be1=0.d0
  alphae1=0.d0
  De1=0.d0
  Be2=0.d0
  alphae2=0.d0
  De2=0.d0
  iprint=0
  rewind(8)
  read(8,nml=basis)
  rmin=rmin!/0.529177d0
  rmax=rmax!/0.529177d0
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
  jmax=0
     do j1=j1min,j1max,j1step
     kmax1=min0(kmax,j1)
!     do ik=0,kmax1
     do ik=-kmax1,kmax1
     kk=iabs(ik)
     kk=mod(kk,3)
     if(kk.EQ.0.AND.kstep.EQ.3)then
     jmax=jmax+1
     elseif(kk.NE.0.AND.kstep.NE.3)then
     jmax=jmax+1
     endif
     enddo
     enddo
  !**********************************************************************************
  write(output,*) ' j1  kap        Ej(cm-1)'
  allocate(k2tab(jmax))
  allocate(ktab(jmax))
  allocate(jtab(jmax))
  allocate(j1tabs(jmax))
  allocate(katabs(jmax))
  !
  !**********************************************************************************
  !
  !   jtab numéro de la voie ----------->  jjtab pour retrouver la voie dans sigma
  !
  !**********************************************************************************
  !
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
     do j1=j1min,j1max,j1step
     kmax1=min0(kmax,j1)
!     do ik=0,kmax1
     do ik=-kmax1,kmax1
     kk=iabs(ik)
     kk=mod(kk,3)
     if(kk.EQ.0.AND.kstep.EQ.3)then
     jmax=jmax+1
        aaa= dfloat(j1*(j1+1))*Be1-dfloat(j1*j1*(j1+1)*(j1+1))*De1 &
             +dfloat(ik*ik)*Bc1
        write(output,'(2I4,2F18.7)') j1,ik, aaa
     flush(output)
        k2tab(jmax)=2.d0*mu*(E-aaa+bbb)/uacmm1
     if(k2tab(jmax)>0.d0) nlevo=nlevo+1
     j1tabs(jmax)=j1
     katabs(jmax)=ik
     elseif(kk.NE.0.AND.kstep.NE.3)then
     jmax=jmax+1
        aaa= dfloat(j1*(j1+1))*Be1-dfloat(j1*j1*(j1+1)*(j1+1))*De1 &
             +dfloat(ik*ik)*Bc1
        write(output,'(2I4,2F18.7)') j1,ik, aaa
     flush(output)
        k2tab(jmax)=2.d0*mu*(E-aaa+bbb)/uacmm1
     if(k2tab(jmax)>0.d0) nlevo=nlevo+1
     j1tabs(jmax)=j1
     katabs(jmax)=ik
     endif
     enddo
     enddo
  !**********************************************************************************
  !
  !       ktab(i) K**2 trié           jtab(i) numéro du niveau i
  !
  !**********************************************************************************
  !
  write(output,*)
  write(output,*) ' n°  j1  kap    k2(cm-1)',nlevo
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
     jj=katabs(i)
     katabs(i)=katabs(m)
     katabs(m)=jj
     write(output,'(3I4,F18.7)') i,j1tabs(i),katabs(i),k2tab(i)*uacmm1/(2.d0*mu)
!     write(output,*)katabs
     flush(output)

  enddo
     write(output,'(3I4,F18.7)') i,j1tabs(i),katabs(i),k2tab(i)*uacmm1/(2.d0*mu)
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
  allocate(ri(0:nmax),hh(0:nmax))
  allocate(fi(0:nmax,0:lamax),fs(0:nmax,0:lamax))
  allocate(m1(0:lamax),l(0:lamax))
  allocate(vvl(0:lamax))

  ri(0)=rmin1!/0.529177d0
  deltar=deltar!/0.529177d0
  do i=1,nmax
     ri(i)=ri(0)+dfloat(i)*deltar
  enddo
  do i=0,nmax-1
     hh(i)=ri(i+1)-ri(i)
  enddo
  fla=.false.
  if(fla)then
     it=0
     do i=0,l1max,2
     do j=0,i
              l(it)=i
              m1(it)=j
              it=it+1
     enddo
     enddo
     do il=0,lamax
        READ(ipot,*) fi(0:nmax,il)
     enddo

  else
     do il=0,lamax
        read(ipot,*) l(il),m1(il)
           READ(ipot,*) fi(0:nmax,il)
     enddo
  endif
!**********************************************************************************
!
!   conversion des Vlambda   (cm-1 ------> ua)
!
!**********************************************************************************
  fi=fi*(2.d0*mu)/uacmm1
!**********************************************************************************
!
!   calcul des dérivés seconde des Vlambda
!
!**********************************************************************************
  do il=0,lamax
     call spline1(ri,fi(:,il),fs(:,il),nmax)
  enddo

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
  do jtot=jtotl,jtotu,jtstep
     !**********************************************************************************
     !
     !............................................. Boucle sur la parité
     !
     !**********************************************************************************
     do parite=0,1
       time_begin = OMP_GET_WTIME()

     !**********************************************************************************
     !
     !................................ calcul du nombre de voies couplées
     !
     !**********************************************************************************
!**********************************************************************************
!
!  détermination du nombre de voies couplées
!
!**********************************************************************************
        nlev=0
        do i=1,jmax   
           j1=j1tabs(i)
           ka=katabs(i)
              lmin=iabs(jtot-j1)
              lmax=jtot+j1
              do il=lmin,lmax
!                 kk=jtot+j1+ka+il+parite
!                 kk=j1+ka+il+parite+mod(ka,3)
!                 kk=j1+ka+il+parite+p6tab(mod(ka,6))
                 kk=j1+ka+il+parite+p6tab(mod(abs(ka),6))
!                 kk=j1+ka+il+parite+mod(mod(abs(ka),3),2)
                 if(ka<0)kk=kk+1
                 if(mod(iabs(kk),2).eq.0) then
                    nlev=nlev+1
                 endif
              enddo
           enddo
        if (nlev/=0) then
           allocate(j1tab(nlev))
           allocate(katab(nlev))
           allocate(jjtab(nlev))
           allocate(ltab(nlev))
           allocate(k2tab(nlev))
           allocate(fl(nlev,nlev,0:lamax))
           allocate(RP(nlev,nlev))
           allocate(ptab(nlev))


!**********************************************************************************
!
!   détermination des voies couplées
!
!**********************************************************************************
           nlev=0
           do i=1,jmax   
           j1=j1tabs(i)
           ka=katabs(i)
                 lmin=iabs(jtot-j1)
                 lmax=jtot+j1
                 do il=lmin,lmax
!                 kk=jtot+j1+ka+il+parite
!                 kk=j1+ka+il+parite+mod(ka,3)
!                 kk=j1+ka+il+parite+p6tab(mod(ka,6))
                 kk=j1+ka+il+parite+p6tab(mod(abs(ka),6))
!                 kk=j1+ka+il+parite+mod(mod(abs(ka),3),2)
                 if(ka<0)kk=kk+1
                 if(mod(iabs(kk),2).eq.0) then
                       nlev=nlev+1
                       j1tab(nlev)=j1
                       katab(nlev)=ka
                       ltab(nlev)=il
                       k2tab(nlev)=ktab(i)
                       jjtab(nlev)=jtab(i)
                       kk=abs(ka)
!                 if(ka<0)kk=kk+1
                       ipar=mod(abs(ka),6)
!                       ipar=0
!                       if(ka<0)ipar=1
                 ptab(nlev)=p6tab(ipar)
                    endif
                 enddo
              enddo
!**********************************************************************************
!
!   calcul de la matrice de couplage  fl
!
!**********************************************************************************
           fl(:,:,:)=0.d0
           call V15(jtot,nlev,fl)
!**********************************************************************************
!
!   propagateur
!
!**********************************************************************************


                 call ijohnson(rmin,rmax,RP,nlev,npas,fl)
!                 call ijohnson_d(rmin,rmax,RP,nlev,npas,fl)
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
              deallocate(ltab,fl,RP,k2tab,j1tab,katab,jjtab,ptab)
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
           Tr=Tr*Tr
           Ti=Ti*Ti
           Tr=(Tr+Ti)*dfloat(2*jtot+1)*dfloat(jtstep)*pi
!**********************************************************************************
!
!   calcul des sections efficaces
!
!**********************************************************************************
           do i=1,nopen
              do j=1,nopen
                 j1=j1tab(i)
                 sigma(jjtab(i),jjtab(j))=sigma(jjtab(i),jjtab(j))+Tr(i,j)/(k2tab(i)*dfloat((2*j1+1)))
              enddo
           enddo
           deallocate(ltab,fl,RP,Tr,Ti,k2tab,j1tab,katab,jjtab,ptab)
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
  write(output,*)'      N°    (  J1,  ka  )          EJ(cm-1)'
do i=1,nlevo
           m=jtab(i)
           j1=j1tabs(i)
           ka=katabs(i)
write(output,2500) i,m,j1,ka, ktab(i)*uacmm1/(2.d0*mu)
enddo
  sigmao(:,:)=sigma(jtab(1:nlevo),jtab(1:nlevo))
  do n=1,nlevo
  do m=1,nlevo
  enddo
  enddo
  call imprime2(iprint,imp,nlevo,sigmao,'sigma     ')
2500   format(2I5,'   (',I4,',',I4,'  )   ',F16.4)
deallocate(jtab,j1tabs,katabs)
end subroutine c0015
