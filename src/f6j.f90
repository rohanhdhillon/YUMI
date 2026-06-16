subroutine f6j(j2,j3,l1,l2,l3,jmin,jmax,r)
  implicit real(kind=8) (a-h,o-z)
integer, parameter :: dp = 8
!integer, parameter :: dp = selected_real_kind(15, 307)
  integer, intent(in) :: j2,j3,l1,l2,l3,jmin,jmax
  integer :: j,jmid
  real (kind=dp) ::  norm, num, den
  real (kind=dp) ::  fj,fj2,fj3,fl1,fl2,fl3
  real (kind=dp), intent(out), dimension(jmin-1:jmax+1) ::  r
  real (kind=dp), dimension(:,:), allocatable ::  a
  integer :: i,jmoins,jplus
  logical :: flagp,flagm
  r=0.d0
  jmid=(jmin+jmax)/2
  allocate(a(jmin-1:jmax+1,3))
  if(jmin.EQ.jmax) then
     r(jmin)=1.d0/dsqrt(dfloat((jmin+jmin+1)*(l1+l1+1)))
     i=mod(j2+j3+l2+l3,2)
     if (i.EQ.1) r(jmin)=-r(jmin)
     return
  endif
  r(:)=0.d0
  a(:,:)=0.d0
  fj2=dfloat(j2)
  fj3=dfloat(j3)
  fl1=dfloat(l1)
  fl2=dfloat(l2)
  fl3=dfloat(l3)
  do j=jmin,jmax
     fj=dfloat(j)
     a(j,1)= dsqrt(((fj+1.d0)*(fj+1.d0))*(((fj*fj)-((fj2-j3)*(fj2-fj3))) &
          *(((fj2+fj3+1.d0)*(fj2+fj3+1.d0))-(fj*fj)) &
          *((fj*fj)-((fl2-fl3)*(fl2-fl3)))*(((fl2+fl3+1.d0)*(fl2+fl3+1.d0))-(fj*fj))))
     a(j,2)=(fj+fj+1.d0)*(fj*(fj+1.d0)*(-fj*(fj+1.d0) &
          +fj2*(fj2+1.d0)+fj3*(fj3+1.d0)-2.d0*fl1*(fl1+1.d0))+fl2*(fl2+1.d0)*(fj*(fj+1.d0)+fj2*(fj2+1.d0) &
          -fj3*(fj3+1.d0))+fl3*(fl3+1.d0)*(fj*(fj+1.d0)-fj2*(fj2+1.d0)+fj3*(fj3+1.d0))) 
     a(j,3)= dsqrt(fj*fj*(((fj+1.d0)*(fj+1.d0)-(fj2-fj3)*(fj2-fj3)) &
          *((fj2+fj3+1.d0)*(fj2+fj3+1.d0)-(fj+1.d0)*(fj+1.d0))*((fj+1.d0)*(fj+1.d0)-(fl2-fl3)*(fl2-fl3)) &
          *((fl2+fl3+1.d0)*(fl2+fl3+1.d0)-(fj+1.d0)*(fj+1.d0))))
  enddo
  do j=jmin,jmax
  enddo
  j=jmax
  flagp=.true.
  jplus=jmid
  jmoins=jmin
  r(jmid)=1.d0
  do while (j.GT.jmid.AND.flagp)
     num=-a(j,1)
     den=a(j,2)+a(j,3)*r(j+1)
     if(dabs(den)-dabs(num).GT.0.d0) then
        r(j)=num/den 
        j=j-1
     else
        jplus=j
        flagp=.false.
        r(j)=1.d0
     endif
  enddo
  do j=jplus,jmoins+1,-1
     r(j-1)=-(a(j,3)*r(j+1)+a(j,2)*r(j))/a(j,1)
  enddo
  norm=r(jmoins)
  if (jmin.EQ.0) then
     r(jmin)=1.d0/dsqrt(dfloat((l2+l2+1)*(j3+j3+1)))
     i=mod(l1+l2+j3,2)
     if (i.EQ.1) r(jmin)=-r(jmin)
     r(jplus)=r(jmin)/norm
     do j=jmin+1,jplus-1
        r(j)=r(j)*r(jplus)
     enddo
     do j=jplus+1,jmax
        r(j)=r(j-1)*r(j)
     enddo
  else
     j=jmin
     flagm=.true.
     jmoins=jplus
     do while (j.LT.jplus.AND.flagm)
        num=-a(j,3)
        den=a(j,2)+a(j,1)*r(j-1)
        if(dabs(den)-dabs(num).GT.0.d0) then
           r(j)=num/den 
           j=j+1
        else
           jmoins=j
           flagm=.false.
           norm=r(j)
           r(j)=1.d0
        endif
     enddo
     r(jplus)=r(jmoins)/norm
     do j=jmin+1,jplus-1
        r(j)=r(j)*r(jplus)
     enddo

     do j=jplus+1,jmax
        r(j)=r(j-1)*r(j)
     enddo
     do j=jmoins-1,jmin,-1
        r(j)=r(j+1)*r(j)
     enddo

  endif
  norm=0.d0
  do j=jmin,jmax
     norm=norm+dfloat(2*j+1)*r(j)*r(j)
  enddo
  norm=dsqrt(norm*dfloat(2*l1+1))
  do j=jmin,jmax
     r(j)=r(j)/norm
  enddo
end subroutine f6j
