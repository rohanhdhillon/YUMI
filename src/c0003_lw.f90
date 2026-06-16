subroutine c0003
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_input
  !$ use omp_lib
  integer :: i, j, kk
  integer :: nlev,lmin,lmax,jmin,jmax,jtot,nopen,nopen1
  real(kind=dp),allocatable, dimension(:,:,:) :: fl
  real(kind=dp),allocatable, dimension(:,:) :: RP,Tr,Ti,sigma,sigmao
  real(kind=dp),allocatable, dimension(:) :: RP0
  real(kind=dp),allocatable, dimension(:,:) :: fitest
  real(kind=dp) :: Be, alphae, De, mu
  real(kind=dp) :: deltar, rmin1
  integer, allocatable, dimension(:) :: j1tabs, j2tabs
  integer :: parite, irest,iprint, imp, output, input8, ipot
  logical :: printw,fla,cuda

  integer :: iligne
  real time_begin,time_end


  namelist /basis/jmin,jmax,jstep,Be1,Be2,alphae1,alphae2,&
       De1,De2,mu,iprint,j1min,j1max,j1step,j2min,j2max,j2step

  !**********************************************************************************
  !
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
  Be1=0.d0
  alphae1=0.d0
  De1=0.d0
  Be2=0.d0
  alphae2=0.d0
  De2=0.d0
  iprint=0
  rewind(8)
  read(8,nml=basis)
  Be1=Be1-alphae1/2.d0
  Be2=Be2-alphae2/2.d0
  mu=mu*FMPRT
  !  mu=mu/BFCT
!**********************************************************************************
!
!   construction de la base
!   détermination de la dimension de la base
!
!**********************************************************************************
  jmax=0
  do j1=j1min,j1max,j1step
     do j2=j2min,j2max,j2step
        jmax=jmax+1
     enddo
  enddo
  write(output,*) ' j1   j2         Ej(cm-1)'
  allocate(k2tab(jmax))
  allocate(ktab(jmax))
  allocate(jtab(jmax))
  allocate(j1tabs(jmax))
  allocate(j2tabs(jmax))
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
  nlevo=0
  jmax=0
  do j1=j1min,j1max,j1step
     do j2=j2min,j2max,j2step
        jmax=jmax+1
        k2tab(jmax)=(2.d0*mu*(E-dfloat(j1*(j1+1))*Be1+dfloat(j1*j1*(j1+1)*(j1+1))*De1))/uacmm1 &
             +(2.d0*mu*(-dfloat(j2*(j2+1))*Be2+dfloat(j2*j2*(j2+1)*(j2+1))*De2))/uacmm1 
        write(output,'(2I4,3F18.7)') j1,j2, dfloat(j1*(j1+1))*Be1-dfloat(j1*j1*(j1+1)*(j1+1))*De1 &
             +dfloat(j2*(j2+1))*Be2-dfloat(j2*j2*(j2+1)*(j2+1))*De2
        if(k2tab(jmax)>0.d0) nlevo=nlevo+1
        j1tabs(jmax)=j1
        j2tabs(jmax)=j2
     enddo
  enddo
  !**********************************************************************************
  !
  !       ktab(i) K**2 trié           jtab(i) numéro du niveau i
  !
  !**********************************************************************************
  !
  write(output,*)
  write(output,*) ' j1   j2         k2(cm-1)'
  do i=1,jmax-1
     m=i
     do j=i+1,jmax
        if(k2tab(j).GT.k2tab(m))m=j
     enddo
     c=k2tab(i)
     k2tab(i)=k2tab(m)
     k2tab(m)=c
     j=jtab(i)
     jtab(i)=jtab(m)
     jtab(m)=j
     j=j1tabs(i)
     j1tabs(i)=j1tabs(m)
     j1tabs(m)=j
     j=j2tabs(i)
     j2tabs(i)=j2tabs(m)
     j2tabs(m)=j
     write(output,'(2I4,F18.7)') j1tabs(i),j2tabs(i),k2tab(i)*219474.63067d0/(2.d0*mu) 

  enddo
  ktab=k2tab
  deallocate(k2tab)
    write(77,'(a,//,a)') ' Results a la Molscat ','   E    Jmin    Jstep     Jmax   ndown   nup   sigma  '

  ! ici est la boucle sur E
!**********************************************************************************
!
!   pas d'intégration
!   (pas utilisé?)
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
  write(*,*) lamax,nmax,rmin1,deltar

  allocate(ri(0:nmax),hh(0:nmax))
  allocate(fi(0:nmax,0:lamax),fs(0:nmax,0:lamax))
  allocate(fitest(0:139,0:lamax))
  allocate(l1(0:lamax),l2(0:lamax),l(0:lamax))
  allocate(vvl(0:lamax))

  ri(0)=rmin1!/0.529177d0
  deltar=deltar!/0.529177d0
  do i=1,nmax
     ri(i)=ri(0)+dfloat(i)*deltar
  enddo
  read(ipot,*,end=999) ri
999 continue
  do i=0,nmax
     write(*,'(i5,g16.6)') i,ri(i)
  end do

  do i=0,nmax-1
     hh(i)=ri(i+1)-ri(i)
  enddo

     do il=0,lamax
        write(*,*) 'reading il = ',il
        read(ipot,*) l1(il),iill,l2(il),l(il)
        write(*,*) l1(il),iill,l2(il),l(il)
        read(ipot,*) fi(0:nmax,il)
        flush(6)
     enddo




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
  flush(output)
  do jtot=jtotl,jtotu,jtstep
     !**********************************************************************************
     !
     !............................................. Boucle sur la parité
     !
     !**********************************************************************************
     j12step=1
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
           j12min=iabs(j1-j2)
           j12max=j1+j2
           do j12=j12min,j12max,j12step
              lmin=iabs(jtot-j12)
              lmax=jtot+j12
              do il=lmin,lmax
                 kk=jtot+1+j1+j2+il+parite
                 if(mod(kk,2).eq.0) then
                    nlev=nlev+1
                 endif
              end do
           enddo
        enddo



        if (nlev/=0) then
           allocate(j1tab(nlev))
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
           j12min=iabs(j1-j2)
           j12max=j1+j2
           do j12=j12min,j12max,j12step
              lmin=iabs(jtot-j12)
              lmax=jtot+j12
              do il=lmin,lmax
                 kk=jtot+1+j1+j2+il+parite
                 if(mod(kk,2).eq.0) then
                    nlev=nlev+1
                       j1tab(nlev)=j1
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
           call V03(jtot,nlev,fl)
!**********************************************************************************
!
!   propagateur
!
!**********************************************************************************
           cuda=.false.
!           if(cuda) then
!                call ijohnson_d(rmin,rmax,RP,nlev,npas,fl)
!           else
                call ijohnson(rmin,rmax,RP,nlev,npas,fl)
!                call rpropagat(rmin,rmax,RP,nlev,npas,fl)
!           endif

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
              deallocate(ltab,fl,RP,k2tab,j1tab,j2tab,j12tab,jjtab,RP0)
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
                 j2=j2tab(i)
                 j1=j1tab(i)
                 sigma(jjtab(i),jjtab(j))=sigma(jjtab(i),jjtab(j))+Tr(i,j)/(k2tab(i)*dfloat((2*j1+1)*(2*j2+1)))
              enddo
           enddo
           deallocate(ltab,fl,RP,Tr,Ti,k2tab,j1tab,j2tab,j12tab,jjtab,RP0)
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
!**********************************************************************************
!
!   voies ouvertes
!
!**********************************************************************************
  write(output,*)
  write(output,*)'      N°    (  J1,  J2)          EJ(cm-1)'
do i=1,nlevo
                 m=jtab(i)
                 j2=j2tabs(i)
                 j1=j1tabs(i)
                 write(output,'(2I5,'   (',I4,',',I4,')   ',F16.4)')
                 & i,m,j1,j2, ktab(m)*uacmm1/(2.d0*mu)
enddo



  sigmao(:,:)=sigma(jtab(1:nlevo),jtab(1:nlevo))
  write(output,*)
  call imprime2(iprint,imp,nlevo,sigmao,'sigma     ')


  do iligne=1,nlevo
      write(


      !  la boucle sur E se termine ici
end subroutine c0003
