subroutine asymm(j1,valp1,vec1,alpha1,epsk1,A,B,C)
! use mod_constantes
!  use mod_pot
!  use mod_base
!  use mod_input
  implicit real*8 (a-h,o-z)
  integer :: i, j, kk
  integer,intent(in) :: j1
  real(kind=8),allocatable, dimension(:,:) :: Hr
  real(kind=8),intent(out), dimension(2*j1+1) :: valp1,alpha1,epsk1
  real(kind=8),intent(out), dimension(2*j1+1,2*j1+1) :: vec1
  real(kind=8),intent(in) :: A,B,C
      REAL(kind=8) , allocatable , dimension(:):: work
      integer(kind=4) :: info
  real(kind=8) :: som

  !**********************************************************************************
  !
  !
  !**********************************************************************************
  !
!  write(*,*)A,B,C
  jp1=2*j1+1
  allocate(Hr(3,jp1))
  allocate(work(3*jp1-2))
  alpha1=1.d0
  epsk1=1.d0
  Hr=0.d0
     i=0
     do  k1=-j1,j1
         i=i+1
         k=iabs(k1)
         k=k1
         Hr(3,i)=0.5d0*(A+B)*dfloat(j1*(j1+1)-k*k)+C*dfloat(k*k)
         if(i<2*j1)Hr(1,i+2)=0.25d0*(A-B)*dsqrt(dfloat((j1*(j1+1)-k*(k+1))*(j1*(j1+1)-(k+1)*(k+2))))
!         if(i>+2)Hr(3,i-2)=0.25d0*(B-C)*dsqrt(dfloat((j1*(j1+1)-k*(k-1))*(j1*(j1+1)-(k-1)*(k-2))))
!         if(i<i-2)Hr(3,i-2)=0.25d0*(B-C)*dsqrt(dfloat((j1*(j1+1)-k*(k-1))*(j1*(j1+1)-(k-1)*(k-2))))
!  write(*,*)j1,k,i,Hr(i,i+2),Hr(i,i-2)
     enddo
!  write(*,*)Hr
      call dsbev('V','U',jp1,2,Hr,3,valp1,vec1,jp1,work,info)
!  write(*,*)
!  write(*,'(15F12.7)')valp1(:)
  do iv=1,jp1
  som=0.d0
  som1=0.d0
  som2=0.d0
  som4=0.d0
  do jv=1,jp1
  som=som+vec1(jv,iv)*vec1(jp1-jv+1,iv)
  som4=som4+vec1(jv,iv)
  ii=mod(jv+j1,2)
  if(ii==0)then
  som1=som1+vec1(jv,iv)*vec1(jv,iv)
  else
  som2=som2+vec1(jv,iv)*vec1(jv,iv)
  endif
  enddo
!  write(*,*)som,som1,som2
!  if(som4<0)vec1(:,iv)=-vec1(:iv)
!  write(*,'(15F12.7)')vec1(:,iv)
  if(som<0.01)epsk1(iv)=-1.d0
  if(dabs(som1)<dabs(som2))alpha1(iv)=-1.d0
  enddo
!  write(*,'(15F12.7)')epsk1(:)
!  write(*,'(15F12.7)')alpha1(:)
  deallocate(Hr,work)
end subroutine asymm
