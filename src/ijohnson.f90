subroutine ijohnson(rmin,rmax,Y,nbv,n,fl,useCholesky)
  !$ use OMP_LIB
  use mod_constantes
  use mod_pot
  use mod_base, only : ltab, jtab, k2tab, lamax, nmax
  implicit real(kind=8) (a-h,o-z)
  integer::INFO,LWORK,LIWORK,ir,lam
  integer , intent(in) :: nbv, n
  real(kind=dp),allocatable,dimension(:,:)::Qa,Qc
  real(kind=dp),allocatable,dimension(:)::r1,r2
  !!!Modified this line of code!!
  real(kind=dp),intent(out),dimension(nbv,nbv)::Y
  real(kind=dp),intent(in),dimension(nbv,nbv,0:lamax)::fl
  integer,allocatable,dimension(:)::isuppz
  real(kind=dp),allocatable,dimension(:)::WORK
  real(kind=dp)::E,k,h,hs3,h2s6,h4,lami
  real(kind=dp), intent(in)::rmin,rmax
  real::result
  logical::iprint
  logical, intent(inout)::useCholesky
  real(kind=dp) :: BIG=1.d20,one=1.0d0,zero=0.d0
  real(kind=8) :: time_begin,time_end


  iprint=.false.
  imp=9
  LWORK=64*nbv
  LIWORK=64*nbv
  h=(rmax-rmin)/dfloat(n)                                      ! Pas d'integration
  hs3=h/3.d0
  h2s6=h*h/6.d0
  h4=4.d0/h
  
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
  
  !$OMP PARALLEL DO PRIVATE (i, j, lam)
  do j=1,nbv
    do lam=0,lamax
      do i=1,j
         Qa(i,j)=Qa(i,j)+fl(i,j,lam)*vvl(lam)
      enddo
     enddo
     Qa(j,j)=Qa(j,j)+dfloat(ltab(j)*(ltab(j)+1))*rm2-k2tab(j)
  enddo
  !$OMP END PARALLEL DO


  !$OMP PARALLEL DO PRIVATE (i , j)
  do j=1,nbv
    Y(j,j)=dsqrt(dabs(Qa(j,j)))
    Qa(j,j)=zero
    do i=1,j
      Qa(i,j)=hs3*Qa(i,j)
    enddo
  enddo
  !$OMP END PARALLEL DO
 
  do ir=1,n,2
    r=r+h                    !*******************  rc
    rm2=1.d0/(r*r)


    !$OMP PARALLEL DO PRIVATE (i , j)
    do j=1,nbv
      do i=1,j
        Qc(i,j)=zero
      enddo
    enddo
    !$OMP END PARALLEL DO

    !Gives vvl via the module mod_pot vvl vector.
    call pot
     

    !$OMP PARALLEL DO PRIVATE (i, j,lam)
    do j=1,nbv
      do lam=0,lamax
        do i=1,j
          !fl comes from v03.f90
          !Here Qc is the potential energy at distance r
          Qc(i,j)=Qc(i,j)+fl(i,j,lam)*vvl(lam)
        enddo
      enddo
    enddo

    !$OMP PARALLEL DO PRIVATE (j) 
    do j=1,nbv
      Qc(j,j)=Qc(j,j)+dfloat(ltab(j)*(ltab(j)+1))*rm2-k2tab(j)
    enddo
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO PRIVATE (i,lami) 
    do i=1,nbv
      lami=dsqrt(dabs(Qc(i,i)))
      !Open channel
      if (Qc(i,i).GT.0.d00) then
        r1(i)=(lami/dtanh(h*lami))
        r2(i)=(lami/dsinh(h*lami))
      !Closed channel
      else
        r1(i)=(lami/dtan(h*lami))
        r2(i)=(lami/dsin(h*lami))
      endif
    enddo
    !$OMP END PARALLEL DO

    !Qa=hs3*w
     
    !$OMP PARALLEL DO PRIVATE (i , j )
    do j=1,nbv
      Qc(j,j)=zero
      do i=1,j
        Qc(i,j)=-h2s6*Qc(i,j)
        Qc(j,i)=Qc(i,j)
      enddo
      Qc(j,j)=Qc(j,j)+one
    enddo
    !$OMP END PARALLEL DO

     
    !Using cholesky decomposition to find the inverse of matrix Qc
    if (UseCholesky) then
      call dpotrf('U',nbv,Qc,nbv,info)
      if (info /=0) then
        write(*,*) "Qc is not positive definite matrix. Cannot use Cholesky Decomposition."
        UseCholesky=.false.
      endif
      
      
      call dpotri('U',nbv,Qc,nbv,info)
      if (info /=0) then
        write(*,*) "Qc is not positive definite matrix. Cannot use Cholesky Decomposition."
        UseCholesky=.false.
      endif
    
    else
    !!Use the following if Cholesky decomposition above does not work.
      if (nbv<500)then
        call dgetrf(nbv,nbv,Qc,nbv,isuppz,info)
        call dgetri(nbv,Qc,nbv,isuppz,work,lwork,info) 
      else
        call matinv(Qc,nbv)
      endif
    endif 
   
    !$OMP PARALLEL DO PRIVATE (i , j)
    do j=1,nbv
      Qc(j,j)=Qc(j,j)-one 
      Qa(j,j)=Qa(j,j)+r1(j)
      do i=1,j
        Qc(i,j)=h4*Qc(i,j)
        Y(i,j)=Y(i,j)+Qa(i,j)
        Qa(i,j)=zero
      enddo
    enddo
    !$OMP END PARALLEL DO

    !
    !
    !!!$OMP SECTIONS
    !
    r=r+h                    !*******************  rb
    !!$OMP SECTION
     
    call pot
     

    !!This must be tested for larger matrices. This is one of the most computationally expensive process.
    !$OMP PARALLEL DO PRIVATE (j, lam) SHARED(fl, vvl)
    do j=1,nbv
      do lam=0,lamax
        do i=1,j-1
          Qa(i,j)=Qa(i,j)+fl(i,j,lam)*vvl(lam)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO PRIVATE (i , j)
    do j=1,nbv
      do i=1,j
        Qa(i,j)=hs3*Qa(i,j)
        Y(j,i)=Y(i,j)
      enddo
    enddo
    !$OMP END PARALLEL DO

     
     !!$OMP SECTION
     call matinv(Y,nbv)
     

     !$OMP PARALLEL DO PRIVATE (i, j)
     do j=1,nbv
       do i=1,nbv
         Y(i,j)=Y(i,j)*r2(i)*r2(j)
       enddo
     enddo
     !$OMP END PARALLEL DO



     !$OMP PARALLEL DO PRIVATE (i , j)
     do j=1,nbv
       Qc(j,j)=Qc(j,j)+r1(j)
       do i=1,j
         Y(i,j)=2*Qc(i,j)-Y(i,j)
         Y(j,i)=Y(i,j)
       enddo
     enddo
     !$OMP END PARALLEL DO

    
     call matinv(Y,nbv)
     


     !$OMP PARALLEL DO PRIVATE (i, j)
     do i=1,nbv
       do j=1,nbv
         Y(j,i)=Y(j,i)*r2(i)*r2(j)
       enddo
     enddo
     !$OMP END PARALLEL DO

     
     !$OMP PARALLEL DO PRIVATE (i , j)
     do j=1,nbv
       do i=1,j
         Y(i,j)=Qa(i,j)-Y(i,j)
         Y(j,i)=Y(i,j)
       enddo
       Y(j,j)=Y(j,j)+r1(j)
     enddo
     !$OMP END PARALLEL DO


     enddo
     !    call imprime2(iprint,imp,nbv,Y,'Y         ')
     deallocate(r1,r2,Qa,Qc)

end subroutine ijohnson
