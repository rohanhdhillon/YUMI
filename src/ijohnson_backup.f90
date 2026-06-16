subroutine ijohnson(rmin,rmax,Y,nbv,n,fl)
     !$ use OMP_LIB
     use mod_constantes
     use mod_pot
     use mod_base, only : ltab, jtab, k2tab, lamax, nmax
     implicit real(kind=8) (a-h,o-z)
     integer :: nbvmax
     integer::INFO,LWORK,LIWORK,ir,lam
     integer , intent(in) :: nbv, n
     real(kind=dp),allocatable,dimension(:,:)::Qa,Qc
     real(kind=dp),allocatable,dimension(:)::r1,r2
     real(kind=dp),intent(out),dimension(nbv,nbv)::Y
     real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
     integer,allocatable,dimension(:)::isuppz
     real(kind=dp),allocatable,dimension(:)::WORK
     real(kind=dp)::E,k,h,hs3,h2s6,h4,lami
     real(kind=dp), intent(in)::rmin,rmax
     real(kind=8)  time_begin,time_end, time_begin1, time_end1, time_begin2, time_end2
     real::result
     real,dimension(2)::tarray
     logical :: iprint
     real(kind=dp) :: BIG=1.d20,one=1.0d0,zero=0.d0
     iprint=.false.
     imp=9
     LWORK=64*nbv
     LIWORK=64*nbv
     h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
     hs3=h/3.d0
     h2s6=h*h/6.d0
     h4=4.d0/h

     nbvmax=500  ! switch matinv/lapack
  
     !  call dtime(tarray, result)                                           ! Allocation de la mémoire
     allocate(r1(nbv),r2(nbv))
     allocate(Qa(nbv,nbv),Qc(nbv,nbv))


     allocate(WORK(LWORK))
     allocate(isuppz(2*nbv))


     r=rmin
     rm2=1.d0/(r*r)
     
     !$OMP PARALLEL DO PRIVATE (i , j )
     do j=1,nbv
       do i=1,nbv
         Y(i,j)=zero
         Qa(i,j)=zero
       enddo
     enddo
     !$OMP END PARALLEL DO
     
     call pot
     
     !OMP PARALLEL DO PRIVATE(j,lam) 
     do lam=0,lamax
       do j=1,nbv
         do i=1,j
           Qa(i,j)=Qa(i,j)+fl(i,j,lam)*vvl(lam)
         enddo
       enddo
     enddo
     !OMP END PARALLEL DO

     !          call imprime2(iprint,imp,nbv,Qa,'Qa        ')

     !$OMP PARALLEL DO PRIVATE (i)
     do i=1,nbv
        Qa(i,i)=Qa(i,i)+dfloat(ltab(i)*(ltab(i)+1))*rm2-k2tab(i)
        Y(i,i)=dsqrt(dabs(Qa(i,i)))
        Qa(i,i)=zero
     enddo
     !$OMP END PARALLEL DO


     !!!$OMP PARALLEL DO PRIVATE (i)  
     !do i=1,nbv
     !   Y(i,i)=dsqrt(dabs(Qa(i,i)))
     !enddo
     !!!$OMP END PARALLEL DO

     !!$OMP PARALLEL DO PRIVATE (i)    
     !!do i=1,nbv
     !!   Qa(i,i)=zero
     !!enddo
     !!$OMP END PARALLEL DO


     !  w=u

     !$OMP PARALLEL DO PRIVATE (i , j)
     do j=1,nbv
        do i=1,j
           Qa(i,j)=hs3*Qa(i,j)
        enddo
     enddo
     !$OMP END PARALLEL DO
     
     do ir=1,n,2
         r=r+h                    !*******************  rc
         rm2=1.d0/(r*r)

         !$OMP PARALLEL DO PRIVATE (i, j)    
         do j=1,nbv
           do i=1,j
             Qc(i,j)=zero
           enddo
         enddo
         !$OMP END PARALLEL DO        
         
         call pot
         
         !$OMP PARALLEL DO PRIVATE(lam,i,j) 
         do j=1,nbv
           do lam=0,lamax 
             do i=1,j
               Qc(i,j)=Qc(i,j)+fl(i,j,lam)*vvl(lam)
             enddo
           enddo
         enddo
         !$OMP END PARALLEL DO
                  
         !$OMP PARALLEL DO PRIVATE (i)
         do i=1,nbv
           Qc(i,i)=Qc(i,i)+dfloat(ltab(i)*(ltab(i)+1))*rm2-k2tab(i)
         enddo
         !$OMP END PARALLEL DO

     
         !$OMP PARALLEL DO PRIVATE (i)
         do i=1,nbv
           lami=dsqrt(dabs(Qc(i,i)))
           if (Qc(i,i).GT.0.d00) then
             r1(i)=(lami/dtanh(h*lami))
             r2(i)=(lami/dsinh(h*lami))
           else
             r1(i)=(lami/dtan(h*lami))
             r2(i)=(lami/dsin(h*lami))
           endif
         enddo
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE (i)
         do i=1,nbv
           Qc(i,i)=zero
         enddo
         !$OMP END PARALLEL DO

         !     Qa=hs3*w
         
         !$OMP PARALLEL DO PRIVATE (i ,j)
         do j=1,nbv
           do i=1,j
             Qc(i,j)=-h2s6*Qc(i,j)
           enddo
         enddo
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE (i ,j)
         do j=1,nbv
           do i=1,j
             Qc(j,i)=Qc(i,j)
           enddo
           Qc(j,j)=Qc(j,j)+one
         enddo
         !$OMP END PARALLEL DO
    

         !call dgetrf(nbv,nbv,Qc,nbv,isuppz,info)
         !call dgetri(nbv,Qc,nbv,isuppz,work,lwork,info)
         if (nbv<nbvmax)then
            call dsytrf('U',nbv,Qc,nbv,isuppz,work,lwork,info)
            call dsytri('U',nbv,Qc,nbv,isuppz,work,info) 
         else         
            call matinv(Qc,nbv)
         endif
  
         !do i=1,nbv
         !!!$OMP PARALLEL DO PRIVATE (j )
         !     do j=i+1,nbv
         !       Qc(j,i)=Qc(i,j)
         !     enddo                                             
         !!!$OMP END PARALLEL DO
         !     enddo

         !$OMP PARALLEL DO PRIVATE (i)
         do i=1,nbv
           Qc(i,i)=Qc(i,i)-one
         enddo
         !$OMP END PARALLEL DO
         
         !$OMP PARALLEL DO PRIVATE (i ,j)
         do j=1,nbv
           do i=1,j
             Qc(i,j)=h4*Qc(i,j)
           enddo
         enddo
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE (i)
         do i=1,nbv
           Qa(i,i)=Qa(i,i)+r1(i)
         enddo
         !$OMP END PARALLEL DO
         
         !$OMP PARALLEL DO PRIVATE (i ,j)
         do j=1,nbv
            do i=1,j
               Y(i,j)=Y(i,j)+Qa(i,j)
            enddo
         enddo
         !$OMP END PARALLEL DO
         !
         !
         !!$OMP SECTIONS
         !
         r=r+h                    !*******************  rb
         !!$OMP SECTION
         !

         !$OMP PARALLEL DO PRIVATE (i ,j)
         do j=1,nbv
           do i=1,j
             Qa(i,j)=zero
           enddo
         enddo 
         !$OMP END PARALLEL DO

         call pot
        
         !$OMP PARALLEL DO PRIVATE (i, j, lam)
         do j=1,nbv
           do lam=0,lamax
             do i=1,j-1
               Qa(i,j)=Qa(i,j)+fl(i,j,lam)*vvl(lam)
             enddo
           enddo
         enddo
         !$OMP END PARALLEL DO

     
         !$OMP PARALLEL DO PRIVATE (i ,j)
         do j=1,nbv
           do i=1,j
             Qa(i,j)=hs3*Qa(i,j)
           enddo
         enddo
         !$OMP END PARALLEL DO

     

         !
         !
         !!$OMP SECTION
         !
         !$OMP PARALLEL DO PRIVATE (i , j )
         do j=1,nbv
           do i=1,j
             Y(j,i)=Y(i,j)
           enddo
         enddo
         !$OMP END PARALLEL DO
        
         !     call dgetrf(nbv,nbv,Y,nbv,isuppz,info)
         !     call dgetri(nbv,Y,nbv,isuppz,work,lwork,info)
   
      
         if (nbv<nbvmax)then
            call dsytrf('U',nbv,Y,nbv,isuppz,work,lwork,info)
            call dsytri('U',nbv,Y,nbv,isuppz,work,info)
         else
            call matinv(Y,nbv)
         endif  
       
        
         !     do i=1,nbv
         !!!$OMP PARALLEL DO PRIVATE (j )
         !     do j=i+1,nbv
         !       Y(j,i)=Y(i,j)
         !     enddo
         !!!$OMP END PARALLEL DO
         !     enddo
         !     call imprime2(iprint,imp,nbv,Y,'Y         ')
         !
         !!$OMP END SECTIONS
         !
        
         time_begin1 = OMP_GET_WTIME()
         
         !!!Suspect the following part of the code is not correct.
         !$OMP PARALLEL DO PRIVATE(i,j)
         do j=1,nbv
           do i=1,nbv
             Y(i,j)=Y(i,j)*r2(i)*r2(j)
           enddo
           Qc(j,j)=Qc(j,j)+r1(j)
         enddo
         !$OMP END PARALLEL DO


         !$OMP PARALLEL DO PRIVATE (i , j)
         do j=1,nbv
           do i=1,j
             Y(i,j)=Qc(i,j)-Y(i,j)
           enddo
         enddo
         !$OMP END PARALLEL DO


     
       !$OMP PARALLEL DO PRIVATE (i , j)
       do j=1,nbv
         do i=1,j
           Y(i,j)=Y(i,j)+Qc(i,j)
         enddo
       enddo
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE (i , j)
       do j=1,nbv
         do i=1,j
           Y(j,i)=Y(i,j)
         enddo
       enddo
       !$OMP END PARALLEL DO
       
       if (nbv<5000)then
          call dsytrf('U',nbv,Y,nbv,isuppz,work,lwork,info)
          call dsytri('U',nbv,Y,nbv,isuppz,work,info)
       else
          call matinv(Y,nbv)
       endif
       
       !$OMP PARALLEL DO PRIVATE (i, j)
       do i=1,nbv
         do j=1,nbv
           Y(j,i)=Y(j,i)*r2(i)
         enddo
       enddo
       !$OMP END PARALLEL DO

       !!!This line of code looks wrong
       !!$OMP PARALLEL DO PRIVATE (i, j )
       !do j=1,nbv
       !  do i=1,nbv
       !    Y(i,j)=Y(i,j)*r2(i)
       !  enddo
       !enddo
       !!$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE (i , j)
       do j=1,nbv
         do i=1,j
           Y(i,j)=Qa(i,j)-Y(i,j)
         enddo
         Y(j,j)=Y(j,j)+r1(j)
       enddo
       !$OMP END PARALLEL DO


       !$OMP PARALLEL DO PRIVATE (i,j)
       do j=1,nbv
         do i=j,nbv
           Y(i,j)=Y(j,i)
         enddo
       enddo
       !$OMP END PARALLEL DO
 
  enddo


 !          call imprime2(iprint,imp,nbv,Y,'Y         ')
deallocate(r1,r2,Qa,Qc)
end subroutine ijohnson
