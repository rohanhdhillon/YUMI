subroutine convergence(sigma,Tr,jmax,nlev,k2tab,jjtab,j1tab,j2tab,nopen,outputFileNum,convergedStatus)
   
  !sigma,Tr,jmax,k2tab,jjtab,j1tab,j2tab,nopen,ConvergedPar1      
  !convergedStatus      
  !$ use OMP_LIB
  
  implicit None
  
  real(kind=8) :: newTerm, oldTerm, tolerance, percentDiff, Difference, maxDiff
  integer,intent(in) :: nopen, jmax, nlev, outputFileNum
  
  real(kind=8),intent(in),dimension(jmax,jmax) :: sigma, Tr
  real(kind=8),intent(in),dimension(nlev) :: k2tab
  integer,intent(in),dimension(nlev) :: jjtab, j1tab, j2tab
  logical,intent(out) :: convergedStatus
  
  !define variables used for the counter
  integer::i,j,j1,j2
 
  tolerance = 1.d-8
  percentDiff = 1.d0
  maxDiff = 1.5d0
   
  
  ! check for the convergence loop
  !!$OMP PARALLEL DO PRIVATE(i, j, j2, j1, oldTerm, newTerm, Difference) SHARED(maxDiff)
  do i=1,nopen
    do j=1,nopen
    
        j2=j2tab(i)
        j1=j1tab(i)
        
        oldTerm = sigma(jjtab(i),jjtab(j))
        newTerm = Tr(i,j)/(k2tab(i)*dfloat((2*j1+1)*(2*j2+1)))
               
        Difference = ABS(newTerm/(OldTerm+tolerance)*100.0)
      
        !store the largest difference
        if (maxDiff<Difference) then
          write(outputFileNum, 1000) 'i', i, "  j ", j, '  oldTerm ', oldTerm, "   newTerm ", newTerm, "  Difference ", Difference
          maxDiff=Difference
        endif
        
          
        if (Difference>percentDiff) then
          exit!!GOTO 100 
        endif
         
      enddo

      !Need two of these to exist both the looks
      if (Difference>percentDiff) then
        exit 
      endif

  enddo
  !100 continue
  !!$OMP END PARALLEL DO
  
  !if everything has difference less than 1 percent of the term. Do this before adding.
  !write(*,*) 'Max Difference is ', maxDiff
  write(outputFileNum,*) 'Max Difference is ', maxDiff
  
  if (maxDiff<percentDiff) then
      convergedStatus = .true.
  endif

  1000 format(A, I7,'  ', A, I7,'  ', A, F14.6,'  ',A,F14.6,'  ',A,F14.6)
  2000 format(A, I7,'  ', A, I7,'  ', A, I7,'  ', A, I7)
end subroutine convergence

