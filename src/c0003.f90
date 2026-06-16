subroutine c0003
   !$ use OMP_LIB
   use mod_constantes
   use mod_pot
   use mod_base
   use mod_input
   implicit none
   integer :: i, j, kk
   integer :: nlev,lmin,lmax,jmin,jmax,jtot,nopen,nclose
   real(kind=dp),allocatable, dimension(:,:,:) :: fl
   real(kind=dp),allocatable, dimension(:,:) :: RP, Tr, Ti, sigma, sigmao
   real(kind=dp),allocatable, dimension(:) :: RP0
   real(kind=dp),allocatable, dimension(:,:) :: fitest
   real(kind=dp) :: Be, alphae, De, mu
   real(kind=dp) :: deltar, rmin1
   integer, allocatable, dimension(:) :: j1tabs, j2tabs
   integer :: parite, irest, iprint, imp, output, input8, ipot, DeleteMe, outputFileNum
   logical :: printw, fla, cuda
   real(kind=8) :: time_begin, time_end
   integer :: j1, j2, nlevo, m, il, iill, j12step, j12max, j12min, j12
   integer :: jstep, j1min, j1max, j1step, j2min, j2max, j2step
   real(kind=4) :: Be2, Be1, c
   real(kind=4) :: sec
   real(kind=dp) :: alphae2, alphae1, De2, De1, drr
  
   !Adding flags to test
   logical :: convergedPar0, convergedPar1

   namelist /basis/jmin,jmax,jstep,Be1,Be2,alphae1,alphae2,&
      De1,De2,mu,iprint,j1min,j1max,j1step,j2min,j2max,j2step

   !These are for saving the file     
   real(kind=4) :: j1ini,j2ini,j1fin,j2fin
   character(len=100) :: filename
   integer :: success 
   logical :: saveData, useCholesky


   !**********************************************************************************
   !                  linear rotor -linear rotor
   !                        v1.5  NJ  26-4-2022
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
   printw=.true.
   saveData=.true.
   useCholesky=.true.
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
   
   
   outputFileNum = 271818
   open(outputFileNum,file="output.dat",status="unknown")
   
   
   !**********************************************************************************
   !
   !   construction of the base
   !   determining the dimension of the base
   !
   !**********************************************************************************

   jmax=0
   do j1=j1min,j1max,j1step
      do j2=j2min,j2max,j2step
        jmax=jmax+1
      enddo
   enddo
  
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
         if(k2tab(jmax)>0.d0) nlevo=nlevo+1
         j1tabs(jmax)=j1
         j2tabs(jmax)=j2
      enddo
   enddo
   !**********************************************************************************
   !
   !       ktab(i) K**2 trié  
   !       jtab(i) numéro du niveau  --- i represents the level number
   !
   !**********************************************************************************
   !
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
   !allocate(fitest(0:139,0:lamax))
   allocate(l1(0:lamax),l2(0:lamax),l(0:lamax))
   allocate(vvl(0:lamax))

   ri(0)=rmin1!/0.529177d0
   deltar=deltar!/0.529177d0
   do i=1,nmax
     ri(i)=ri(0)+dfloat(i)*deltar
   enddo

   read(ipot,*) ri
   do i=0,nmax-1
     hh(i)=ri(i+1)-ri(i)
   enddo

   do il=0,lamax
     read(ipot,*) l1(il),iill,l2(il),l(il)
     read(ipot,*) fi(0:nmax,il)
   enddo


   !***************************************************************
   !       conversion des (of) Vlambda   (cm-1 ------> ua)        !
   !***************************************************************
   fi=fi*(2.d0*mu)/uacmm1

   !***************************************************************
   !   calcul des dérivés seconde des Vlambda                     !
   !   calculate the second derivative of Vlambda                 !
   !***************************************************************
   do il=0,lamax
     call spline1(ri,fi(:,il),fs(:,il),nmax)
   enddo

   allocate(sigma(jmax,jmax))
   allocate(sigmao(nlevo,nlevo))
   sigma=0.d0
   
   
   write(output,*)'                      voies          voies ' 
   write(output,*)'      J  par         totales        ouvertes       Temps CPU'
   flush(output)

   
   write(outputFileNum,*) "The jtotl", jtotl, "  jtotu is: ", jtotu, "  jtstep is: ", jtstep

   
   
   !****************************************
   !          Boucle sur J total           !
   !          Loop over the JTot           !
   !****************************************
   
   !Make sure both the converged parities values of the sigma value have converged

   do jtot=jtotl,jtotu,jtstep   
      j12step=1
      convergedPar0=.false.
      convergedPar1=.false.
      do parite=0,1  ! Boucle sur la parité-- Looping over the parity
         time_begin = OMP_GET_WTIME()

         !******************************************************************!
         !        détermination du nombre de voies couplée                  !
         !        Determine the number of coupled channels                  !
         !******************************************************************!
         
         nlev=0
         do i=1,jmax
            j1=j1tabs(i)
            j2=j2tabs(i)
            j12min=iabs(j1-j2)
            j12max=j1+j2
            
            do j12=j12min,j12max,j12step
               
               !il is the orbital angular momentum
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

            !****************************************************************************
            !                  détermination des voies couplées                         !   
            !                 determination of the coupled channels                     !                                                                             !               
            !****************************************************************************
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
            !*********************************************************!
            !       calcul de la matrice de couplage  fl              !
            !         Calculate the coupling matrix fl                !
            !**********************************************************
            
            time_begin = OMP_GET_WTIME()
            call V03(jtot,nlev,fl)
            time_end = OMP_GET_WTIME()

           
            write(outputFileNum, '(A,I9.7F9.2)') "Time taken for V03 is: ", nlev,time_end-time_begin  
            flush(outputFileNum)  
            
            !*********************************************
            !         propagateuri / propagator          !
            !*********************************************
            
            time_begin = OMP_GET_WTIME()
            call ijohnson(rmin,rmax,RP,nlev,npas,fl, useCholesky)
            time_end = OMP_GET_WTIME()

            write(outputFileNum,'(A,I9.7F9.2)') "Time taken for ijohnson: ",nlev, time_end-time_begin
              
            !$OMP PARALLEL DO PRIVATE(i,j)
            do i=1,nlev
               do j=i,nlev
                  RP(i,j)=(RP(i,j)+RP(j,i))/2.d0
                  RP(j,i)=RP(i,j)
               enddo
            enddo
            !$OMP END PARALLEL DO 
                      
            !**********************************************************************************!  
            !             détermination du nombre de voies couplées ouvertes                   ! 
            !          determination of the number of the open coupled channels                !              
            !**********************************************************************************!
           
            nopen = 0
            !!!$OMP PARALLEL DO PRIVATE(j) REDUCTION(+:nopen)
            do j=1,nlev
               if(k2tab(j)>0.d0) then
                  nopen=nopen+1
               endif
            enddo
            !!!$OMP END PARALLEL DO 
            
            if(nopen.EQ.0) then
               deallocate(ltab,fl,RP,k2tab,j1tab,j2tab,j12tab,jjtab,RP0)
               cycle
            endif

            !*******************************************************************************!
            !              calcul de la matrice    T  = Tr  +  i Ti                         !
            !              Calculate the matrix    T =  Tr  +  i Ti                         !
            !*******************************************************************************!
           
            !nopen1=nopen
            allocate(Tr(nopen,nopen), Ti(nopen,nopen))
            call tmat(RP,nlev,nopen,ltab,k2tab,rmax,Tr,Ti,iprint,imp)
           
            !!!Save both Tr, Ti and all details: partity, Total J, Kinetic Energy and Quantum Numbers, nopen1, j1tab, j2tab, and k2tab  
           
            if (saveData) then
               
               !!!Saving Tr matrix
               write(filename, '(A,I3.3,A,I3.3,A,I5.5,A)') 'Tr_JTot_', jtot,"_Par_",parite,"_Size_",nopen,".dat"
               open(unit=31415, file=filename, status='replace' , form="unformatted", action='write', iostat=success)
               write(31415) Tr
               close(31415)


               !!!Saving Ti matrix
               write(filename, '(A,I3.3,A,I3.3,A,I5.5,A)') 'Ti_JTot_', jtot,"_Par_",parite,"_Size_",nopen,".dat"
               open(unit=31416, file=filename, status='replace', form="unformatted", action='write', iostat=success)
               write(31416) Ti   
               close(31416)
              
            endif          

            !OMP PARALLEL DO PRIVATE(i,j)
            do i=1,nopen
               do j=1,nopen
                  Tr(i,j)=(Tr(i,j)**2+Ti(i,j)**2)*dfloat(2*jtot+1)*dfloat(jtstep)*pi
               enddo
            enddo
            !OMP END PARALLEL DO

            !***********************************************************!
            !              calcul des sections efficaces                !
            !          Calculation of the effective cross-sections      !
            !***********************************************************!         
            if (parite==0) then
               !call convergence(sigma,Tr,jmax,nlev,k2tab,jjtab,j1tab,j2tab,nopen,ConvergedPar0)
               if (jtot>5) then
                  !!FIND ME what is the difference between nlevo and nopen
                  call convergence(sigma,Tr,jmax,nlev,k2tab,jjtab,j1tab,j2tab,nopen, outputFileNum, ConvergedPar0)
                  write(outputFileNum,*) "jtot:", jtot, "Par: ", parite,"    ConvergedStatus", ConvergedPar0
                  flush(outputFileNum)
               end if
            else
               if (jtot>5) then
                  call convergence(sigma,Tr,jmax,nlev,k2tab,jjtab,j1tab,j2tab,nopen, outputFileNum, ConvergedPar1)
                  write(outputFileNum,*) "jtot:", jtot, "Par: ", parite,"    ConvergedStatus", ConvergedPar1
                  flush(outputFileNum)
               end if
            endif
            write(outputFileNum,*) "   "
            write(outputFileNum,*) "   "


            !OMP PARALLEL DO PRIVATE(i,j)
            do i=1,nopen
               do j=1,nopen
                  j2=j2tab(i)
                  j1=j1tab(i)
                  sigma(jjtab(i),jjtab(j))=sigma(jjtab(i),jjtab(j))+Tr(i,j)/(k2tab(i)*dfloat((2*j1+1)*(2*j2+1)))
               enddo
            enddo
            !OMP END PARALLEL DO
         
            !Deallocating the matrix
            deallocate(ltab,fl,RP,Tr,Ti,k2tab,j1tab,j2tab,j12tab,jjtab,RP0)

            !exit if both of them have converge
            if (ConvergedPar0.and.ConvergedPar1) then
               write(outputFileNum,*) "Both the parities have converged"
               GOTO 300
            endif

         else
            nopen=0
         endif

         time_end = OMP_GET_WTIME()

         write(output,'(I8,I3,2I15,F19.3)') jtot,parite,nlev,nopen,time_end-time_begin
         flush(output)
      enddo   !---> fin parite/end loop for parity
   enddo      !---> fin J total/enddo for J total 
   

   300 continue  
   !*************************************************
   !       conversion Angsetrum**2 <--- ua**2       !
   !    convert the units Angstrom**2 <--- ua**2    !
   !*************************************************
   sigma=sigma*A02A2

   
   if (saveData) then      
      !!!Saving Smatrix
      write(filename, '(A,I5.5,A)') "SigmaMatrix_Energy_Size_",nopen,".dat" 
      open(unit=31517, file=filename, form="unformatted", action='write', iostat=success)
      write(31517) sigma
      
      close(31517)
   endif

   !********************************************
   !   voies ouvertes/open channels            !
   !********************************************
   write(output,*)
   write(output,*)'      N°    (  J1,  J2)          EJ(cm-1)'
   do i=1,nlevo
      m=jtab(i)
      j2=j2tabs(i)
      j1=j1tabs(i)
      write(output,2500) i,m,j1,j2, ktab(m)*uacmm1/(2.d0*mu)
   enddo
   sigmao(:,:)=sigma(jtab(1:nlevo),jtab(1:nlevo))

   if (saveData) then
      !!!Saving Smatrix
      write(filename, '(A,I5.5,A)') "SMatrix_Energy_Size_",nlevo,".dat"
      open(unit=31418, file=filename, form="unformatted", action='write', iostat=success)
      write(31418) sigmao
      close(31418)
   endif

   write(output,*)
   call imprime2(iprint,imp,nlevo,sigmao,'sigma     ')
   2500   format(2I5,'   (',I4,',',I4,')   ',F16.4)
end subroutine c0003
