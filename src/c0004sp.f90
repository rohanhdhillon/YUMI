subroutine c0004sp
 use mpi
 use mod_constantes
  use mod_pot
  use mod_base
  use mod_input
!$ use omp_lib        
  integer :: i, j, kk,isym
  integer :: nlev,lmin,lmax,jmin,jmax,jtot,nopen,nopen1
!  real(kind=dp),allocatable, dimension(:) :: flsp
!  integer,allocatable, dimension(:) :: indexl,indexj
  real(kind=dp),allocatable, dimension(:,:) :: vec1
  real(kind=dp),allocatable, dimension(:,:) :: RP,Tr,Ti,sigma,sigmao,sigmas,sigmat
  integer , allocatable, dimension(:) :: j1tabs, katabs, j2tabs
  integer , allocatable, dimension(:,:) :: jvaleurs
  real(kind=dp),allocatable, dimension(:) :: RP0
  real(kind=dp) :: Be, alphae, De, mu,A,B,C,Emax
  real(kind=dp) :: deltar, rmin1
  integer :: parite, irest,iprint, imp, output, input8, ipot,icase
  logical :: printw,fla
  real time_begin,time_end
  integer :: nb_procs,rang,longueur_tranche,code
  integer, dimension(MPI_STATUS_SIZE) :: statut
  namelist /basis/jmin,jmax,jstep,Be1,Be2,alphae1,alphae2,&
       De1,De2,mu,iprint,j1min,j1max,j1step,j2min,j2max,k1max,j2step,Bc1,&
       isym,A,B,C,Emax

!**********************************************************************************
!               Toupie symétrique -  Rotateur Linéaire   
!
!      Phillips et al. J. Chem. Phys. 102, 6024 (1995); doi: 10.1063/1.469337 
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
  is11=11
  itmat=15
  printw=.false.
  jtot=3
  j1min=0
  j1max=3
  j1step=1
  j2min=0
  j2max=3
  j2step=1
  k1max=j1max
  isym=1
  Be1=0.d0
  alphae1=0.d0
  De1=0.d0
  Be2=0.d0
  alphae2=0.d0
  De2=0.d0
  A=1.d0
  B=1.d0
  C=1.d0
  Emax=10000.d0
  iprint=0
  icase=1
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
     valp=0.d0
     vec=0.d0
     do j1=j1min,j1max,j1step
        jp1=j1+j1+1
        call asymm(j1,valp(j1,1:jp1),vec(j1,1:jp1,1:jp1),alpha(j1,1:jp1),epsk(j1,1:jp1),A,B,C)
     enddo
  jmax=0
  do j2=j2min,j2max,j2step
     do j1=j1min,j1max,j1step
        jp1=j1+j1+1
     do ik=1,jp1
!     ipar=mod((ik+j1+1),2)
!     ipar=mod(ik/2+(jp1-ik+1)/2,2)
if (icase==1)then
     ipar=mod((jp1-ik+1)/2,2)
     ipar=mod((ik/2),2)
!     if(ipar==isym1)jmax=jmax+1
else
     ipar=mod((ik+j1+1),2)
endif
!     ipar=mod((int(epsk(j1,ik))+1),2)
!     ipar=1
!     if(int(epsk(j1,ik))==+1)ipar=0
!     ipar=mod(ipar+j1,2)
!     if(ipar==isym)jmax=jmax+1
     if(ipar==isym)then
        aaa= valp(j1,ik) &
             +dfloat(j2*(j2+1))*Be2-dfloat(j2*j2*(j2+1)*(j2+1))*De2
     if(aaa<Emax)     jmax=jmax+1
     endif
     enddo
     enddo
  enddo
  !**********************************************************************************
  write(output,*) ' j2  j1  ka   kc  tau  par      Ej(cm-1)',jmax,' voies totales'
  allocate(k2tab(jmax))
  allocate(ktab(jmax))
  allocate(jtab(jmax))
  allocate(j1tabs(jmax))
  allocate(j2tabs(jmax))
  allocate(katabs(jmax))
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
     do j1=j1min,j1max,j1step
        jp1=j1+j1+1
     do ik=1,jp1
!     ipar=mod(ik/2+(jp1-ik+1)/2,2)
!!!!!!!     ipar=mod((jp1-ik+1)/2,2)
!!!!!!!     ipar=mod((ik/2),2)
!     ipar=1
!     if(int(epsk(j1,ik))==+1)ipar=0
!     ipar=mod(ipar+j1,2)
if (icase==1)then
     ipar=mod((jp1-ik+1)/2,2)
     ipar=mod((ik/2),2)
!     if(ipar==isym1)jmax=jmax+1
else
     ipar=mod((ik+j1+1),2)
endif
     if(ipar==isym)then
        aaa= valp(j1,ik) &
             +dfloat(j2*(j2+1))*Be2-dfloat(j2*j2*(j2+1)*(j2+1))*De2
     if(aaa<Emax)     jmax=jmax+1
!        if(jmax==1)bbb=aaa
        if(aaa<Emax)then
        write(output,'(6I4,2F18.7)') j2,j1,(ik)/2,(jp1-ik+1)/2,ik,int(epsk(j1,ik)), aaa
     flush(output)
        k2tab(jmax)=2.d0*mu*(E-aaa+bbb)/uacmm1
     if(k2tab(jmax)>0.d0) nlevo=nlevo+1
     j1tabs(jmax)=j1
     j2tabs(jmax)=j2
     katabs(jmax)=ik
     endif
     endif
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
  write(output,*) ' j2   j1  ka  kc        k2(cm-1)',nlevo,' voies ouvertes'
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
     write(output,'(4I4,F18.7)') j2tabs(i),j1tabs(i),katabs(i)/2,(2*j1tabs(i)-katabs(i)+2)/2,k2tab(i)*uacmm1/(2.d0*mu) 
!     write(output,*)katabs
     flush(output)

  enddo
  i=jmax
     write(output,'(4I4,F18.7)') j2tabs(i),j1tabs(i),katabs(i)/2,(2*j1tabs(i)-katabs(i)+2)/2,k2tab(i)*uacmm1/(2.d0*mu) 
  ktab=k2tab
  deallocate(k2tab)

!**********************************************************************************
!
!   pas d'intégration
!
!**********************************************************************************
!  drr=(pi/dsqrt(E*(2.d0*mu)/uacmm1))/dble(npas)
  drr=(pi/dsqrt(ktab(1)))/dble(npas*2)
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
  allocate(l1(0:lamax),m1(0:lamax),l2(0:lamax),l(0:lamax))
  allocate(vvl(0:lamax))

    irf=1
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
        read(ipot,*) l1(il),m1(il),l2(il),l(il)
!        write(output,*) l1(il),m1(il),l2(il),l(il),lamax
           READ(ipot,*) fi(0:nmax,il)
!           write(OUTPUT,*) fi(irest*100:nmax,il)
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
  allocate(sigmas(jmax,jmax))
  sigma=0.d0
  !**********************************************************************************
  !
  !............................................. Boucle sur J total
  !
  !**********************************************************************************
  write(output,*)'                      voies          voies ' 
  write(output,*)'      J  par         totales        ouvertes       Temps CPU'
  flush(output)
     j12step=1
     i=0
     



call MPI_INIT (code)
call MPI_COMM_SIZE ( MPI_COMM_WORLD ,nb_procs,code)
call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)
nb_tranche=(jtotu-jtotl+1)/nb_procs
if (rang == 0) then
allocate(jvaleurs(nb_procs,nb_tranche))
allocate(sigmat(nlevo,nlevo))
  sigmat=0.d0
do j=1,nb_tranche
do i=1,nb_procs
jvaleurs(i,j)=i-1+(j-1)*nb_procs
enddo
enddo
end if
do j_tranche=1,nb_tranche
call MPI_SCATTER (jvaleurs(:,j_tranche),1, MPI_INTEGER ,jtot,1, &
MPI_REAL ,0, MPI_COMM_WORLD ,code)





!  do jtot=jtotl,jtotu,jtstep
     !**********************************************************************************
     !
     !............................................. Boucle sur la parité
     !
     !**********************************************************************************
     do parite=0,1
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
           j12min=iabs(j1-j2)
           j12max=j1+j2
           do j12=j12min,j12max,j12step
              lmin=iabs(jtot-j12)
              lmax=jtot+j12
              do il=lmin,lmax
                 kk=jtot+j1+j2+il+parite
                 if(epsk(j1,ka)>0)kk=kk+1
                 kk=iabs(kk)
!                 kk=0
                    if(mod(kk,2).eq.0) then
                    nlev=nlev+1
                    endif
              end do
           enddo
        enddo
        if (nlev/=0) then
           allocate(j1tab(nlev))
           allocate(katab(nlev))
           allocate(j2tab(nlev))
           allocate(j12tab(nlev))
           allocate(jjtab(nlev))
           allocate(ltab(nlev))
           allocate(k2tab(nlev))
!  allocate(flsp(0),indexl(0),indexj(0))
           allocate(RP(nlev,nlev))

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
              j12min=iabs(j1-j2)
              j12max=j1+j2
              do j12=j12min,j12max,j12step
                 lmin=iabs(jtot-j12)
                 lmax=jtot+j12
                 do il=lmin,lmax
                 kk=jtot+j1+j2+il+parite
                 if(epsk(j1,ka)>0)kk=kk+1
                 kk=iabs(kk)
!                 kk=0
                    if(mod(kk,2).eq.0) then
                       nlev=nlev+1
                       j1tab(nlev)=j1
                       katab(nlev)=ka
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
           call V04sp(jtot,nlev,nsp)


!**********************************************************************************
!
!   propagateur
!
!**********************************************************************************
                call ijohnsonsp(rmin,rmax,RP,nlev,npas,nsp)
!           do i=1,nlev
!              do j=i,nlev
!                 RP(i,j)=(RP(i,j)+RP(j,i))/2.d0
!                 RP(j,i)=RP(i,j)
!              enddo
!           enddo
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
              deallocate(ltab,flsp,indexl,indexj,RP,k2tab,j1tab,katab,j2tab,j12tab,jjtab)
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
           deallocate(ltab,flsp,indexl,indexj,RP,Tr,Ti,k2tab,j1tab,katab,j2tab,j12tab,jjtab)
        else
           nopen=0
        endif
  time_end = OMP_GET_WTIME() 
!!!  write(output,'(I8,I3,2I15,F19.3)') jtot,parite,nlev,nopen,time_end-time_begin
!call MPI_COMM_RANK ( MPI_COMM_WORLD ,rang,code)
  print *, jtot,parite,nlev,nopen,time_end-time_begin
  flush(output)
     enddo
     !**********************************************************************************
     !
     !..........................................     fin parité
     !
     !**********************************************************************************
!  enddo
!print *,nlevo,jmax,rang
!print *,sigma
!ntot=jmax*jmax
!call MPI_REDUCE (sigma,sigmas,ntot, MPI_REAL , MPI_SUM ,0, MPI_COMM_WORLD ,code)

!print *,nlevo
!  sigmao(:,:)=sigmas(jtab(1:nlevo),jtab(1:nlevo))
!yyprint *,sigmao
!  call imprime2(iprint,imp,nlevo,sigmao,'sigma     ')
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
!  sigmas=sigmas*A02A2
!!!  write(output,*)
!  write(output,*)'      N°    (  J1,  ka,  J2)          EJ(cm-1)'
!!!  write(output,*) ' N°  j2   j1  ka  kc        k2(cm-1)'
!do i=1,nlevo
!           m=jtab(i)
!           j1=j1tabs(i)
!           j2=j2tabs(i)
!           ka=katabs(i)
!!!     write(output,'(5I4,F18.7)') i,j2,j1,ka/2,(2*j1-ka+2)/2,ktab(i)*uacmm1/(2.d0*mu) 
!write(output,2500) i,m,j1,ka,j2, ktab(i)*uacmm1/(2.d0*mu)
!enddo
!call MPI_BARRIER ( MPI_COMM_WORLD ,code)
enddo
ntot=jmax*jmax
call MPI_REDUCE (sigma,sigmas,ntot, MPI_DOUBLE_PRECISION , MPI_SUM ,0, MPI_COMM_WORLD ,code)
!print *,nlevo,jmax,rang
!print *,sigma(1,1)*A02A2
if (rang == 0) then
  sigmas=sigmas*A02A2
  sigmao(:,:)=sigmas(jtab(1:nlevo),jtab(1:nlevo))
!  write(output,*)
!  write(is11)sigmao
!  call imprime2(iprint,imp,nlevo,sigmao,'sigma     ')
print *,nlevo,jmax,rang
sigmat=sigmat+sigmao
do j=1,nlevo
print *,sigmat(:,j)
enddo
endif
2500   format(2I5,'   (',I4,',',I4,',',I4,')   ',F16.4)
call MPI_FINALIZE (code)
end subroutine c0004sp
