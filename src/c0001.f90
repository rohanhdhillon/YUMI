subroutine c0001
  !integer, parameter :: dp = selected_real_kind(15, 307)
  use mod_constantes
  use mod_pot
  use mod_base
  use mod_input
!$ use omp_lib
  !        use cudafor
  implicit real(kind=dp) (a-h,o-z)
  integer :: i, j, kk
  integer :: nlev,lmin,lmax,jmin,jmax,jtot,nopen,nopen1
  real, dimension(2) :: tarray
  real :: result
  real(kind=dp),allocatable, dimension(:,:) :: RP,Tr,Ti,sigma,sigmao
  real(kind=dp),allocatable, dimension(:) :: bid
  real(kind=dp) :: Be, alphae, De, mu
  real(kind=dp) :: deltar,  rmin1
  real(kind=dp),allocatable,dimension(:,:,:)::fl
  integer :: parite, irest,iprint, imp, output, input8, ipot
  logical :: printw,cuda
  real time_begin,time_end
  namelist /basis/jmin,jmax,jstep,Be,alphae,&
       De,mu,iprint,cuda,S2
  real :: time, random
  integer :: istat


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
  !**********************************************************************************
  !
  !   construction de la base
  !   détermination de la dimension de la base
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   pas d'intégration
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   lecture des Vlambda
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   conversion des Vlambda   (cm-1 ------> ua)
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   calcul des dérivés seconde des Vlambda
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !  détermination du nombre de voies couplées
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   détermination des voies couplées
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   calcul de la matrice de couplage  fl
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   propagateur
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !  détermination du nombre de voies couplées ouvertes
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   calcul de la matrice    T  = Tr  +  i Ti
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   calcul des sections efficaces
  !
  !**********************************************************************************
  !**********************************************************************************
  !
  !   conversion Angsetrum**2 <--- ua**2
  !
  !**********************************************************************************
  jtot=3
  jmin=0
  jmax=3
  Be=0.d0
  alphae=0.d0
  De=0.d0
  iprint=0
  cuda=.false.
  S2=0.d0
  rewind(8)
  read(8,nml=basis)
  !**********************************************************************************
  !
  !   construction de la base
  !
  !**********************************************************************************
  Be=Be-alphae/2.d0
  mu=mu*fmprt
  write(output,*) ' j        Ej'
  do j=0,jmax
     write(output,'(I4,F18.7)') j, dfloat(j*(j+1))*Be-dfloat(j*j*(j+1)*(j+1))*De
  enddo
  nlevo=0
  do j=0,jmax
     ak2=(E-dfloat(j*(j+1))*Be+dfloat(j*j*(j+1)*(j+1))*De)
     if(ak2>0.d0) nlevo=nlevo+1
  enddo
  nlevo=nlevo-1
  allocate(sigmao(0:nlevo,0:nlevo))

  !**********************************************************************************
  !
  !   pas d'intégration
  !
  !**********************************************************************************
  drr=pi/dsqrt(E*(2.d0*mu)/219474.63067d0)/npas
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

  read(ipot,*) lamax,nmax,rmin1,deltar
  allocate(ri(0:nmax),hh(0:nmax))
  allocate(fi(0:nmax,0:lamax),fs(0:nmax,0:lamax))
  ri(0)=rmin1!/0.529177d0
  deltar=deltar!/0.529177d0
  do i=1,nmax
     ri(i)=ri(0)+dfloat(i)*deltar
  enddo
  !  rmax=ri(nmax)
  do i=0,nmax-1
     hh(i)=ri(i+1)-ri(i)
  enddo

  do il=0,lamax
     read(ipot,*) i
     !     write(*,*)i
     irest=((nmax+1)/550)
     !     write(*,*)i,irest,nmax
     do j=1,irest
        READ(ipot,*) fi((j-1)*550:j*550-1,il)
        !        write(*,*) fi((j-1)*150:j*150-1,il)
     enddo
     if (irest*550 .NE. nmax+1) then
        READ(ipot,*) fi(irest*550:nmax,il)
        !        write(*,*) fi(irest*550:nmax,il)
     endif
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

  allocate(sigma(0:jmax,0:jmax))
  allocate(vvl(0:lamax))
  sigma=0.d00
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
        !  write(2,*)'jtot,j,l,kk,lmin,lmax'
        !**********************************************************************************
        !
        !  détermination du nombre de voies couplées
        !
        !**********************************************************************************
        nlev=0
        do j=jmin,jmax,jstep
           lmin=iabs(jtot-j)
           lmax=jtot+j
           do il=lmin,lmax
              kk=j+il+parite
              if((((kk/2)*2)-kk).eq.0) then
                 nlev=nlev+1
              endif
           end do
        enddo
        if (nlev/=0) then
           allocate(jtab(nlev))
           allocate(ltab(nlev))
           allocate(k2tab(nlev))
           allocate(RP(nlev,nlev))
           allocate(fl(nlev,nlev,0:lamax))

           !**********************************************************************************
           !
           !   détermination des voies couplées
           !
           !**********************************************************************************


           nlev=0
           do j=jmin,jmax,jstep
              lmin=iabs(jtot-j)
              lmax=jtot+j
              do il=lmin,lmax
                 kk=j+il+parite
                 if((((kk/2)*2)-kk).eq.0) then
                    nlev=nlev+1
                    jtab(nlev)=j
                    ltab(nlev)=il
                 endif
              end do
           enddo
           do i=1,nlev
              j=jtab(i)
              k2tab(i)=(2.d0*mu*(E-dfloat(j*(j+1))*Be+dfloat(j*j*(j+1)*(j+1))*De))/219474.63067d0 
           enddo
           !**********************************************************************************
           !
           !   calcul de la matrice de couplage  fl
           !
           !**********************************************************************************
           call Percival(jtot,nlev,fl)

           !**********************************************************************************
           !
           !   propagateur
           !
           !**********************************************************************************


           flush(output)
           !           if(cuda) then
           !                call ijohnson_d(rmin,rmax,RP,nlev,npas,fl)
           !           else
           call ijohnson(rmin,rmax,RP,nlev,npas,fl)
           !                call rpropagat(rmin,rmax,RP,nlev,npas,fl)
           !           endif
           !             call imprime2(iprint,imp,nlev,RP,'RP        ')

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
              deallocate(ltab,fl,RP,k2tab,jtab)
              cycle
           endif
           !**********************************************************************************
           !
           !   calcul de la matrice    T  = Tr  +  i Ti
           !
           !**********************************************************************************
           !nopen1=nlev
           nopen1=nopen
           !nopen1=nlev
           allocate(Tr(nopen1,nopen1),Ti(nopen1,nopen1))
           allocate(bid((nopen1*(nopen1+1))/2))
           call tmat(RP,nlev,nopen1,ltab,k2tab,rmax,Tr,Ti,iprint,imp)
           Tr=Tr*Tr
           Ti=Ti*Ti
           Tr=(Tr+Ti)*dfloat(2*jtot+1)*dfloat(jtstep)*pi
           !  call imprime2(iprint,imp,nopen1,Tr,'T2        ')
           !**********************************************************************************
           !
           !   calcul des sections efficaces
           !
           !**********************************************************************************
           do i=1,nopen
              do j=1,nopen
                 sigma(jtab(i),jtab(j))=sigma(jtab(i),jtab(j))+Tr(i,j)/(k2tab(i)*dfloat(2*jtab(i)+1))
              enddo
           enddo
           !  call imprime2(iprint,imp,jmax+1,sigma,'sigma     ')
           deallocate(jtab,ltab,fl,RP,Tr,Ti,k2tab,bid)
        else
           nopen=0
        endif
        time_end = OMP_GET_WTIME() 
        write(output,'(I8,I3,2I15,F19.3)') jtot,parite,nlev,nopen,time_end-time_begin
        call flush(output)
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
  sigmao=sigma(0:nlevo,0:nlevo)*A02A2
  write(output,*)  
  write(output,*)'*********************'  
  write(output,*)  
  call imprime2(iprint,imp,nlevo+1,sigmao,'sigma     ')
end subroutine c0001
