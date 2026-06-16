function f3j0(jj1,jj2,jj3)
  implicit real(kind=8) (a-h,o-z)
integer, parameter :: dp = 8
!integer, parameter :: dp = selected_real_kind(15, 307)
  integer,intent(in) :: jj1,jj2,jj3
  integer :: j,j1,j2,j3
  real (kind=dp) ::  f
  real (kind=dp) ::  f3j0
  j1=jj1
  j2=jj2
  j3=jj3
  j=((j1+j2+j3)/2)*2-(j1+j2+j3)
  if(j.ne.0)then
     f3j0=0.d0
     return
  endif
  if(j3.lt.j2)then
     j=j2
     j2=j3
     j3=j
  endif
  if(j2.lt.j1)then
     j=j1
     j1=j2
     j2=j
  endif
  if(j3.lt.j2)then
     j=j2
     j2=j3
     j3=j
  endif
  j=(j1+j2+j3)/2
  f=dfloat(2*j+1)
  do i=1,j1
     f=f*dfloat((2*j-2*i+1)*2*i)/dfloat((2*j-2*i+2)*(2*i-1))
  enddo
  do i=1,j1+j2-j
     f=f*dfloat((2*j1-2*i+1)*2*i)/dfloat((2*j1-2*i+2)*(2*i-1))
  enddo
  f=1.d0/dsqrt(f)
  j=j-(j/2)*2
  if(j.ne.0) f=-f
  f3j0=f
  return
end function f3j0
