      function dlmn(l,m,n,theta)
               implicit real(kind=8) (a-h,o-z) 
               integer, parameter :: dp = 8
               real(kind=dp) :: dlmn,sinb,cosb
               real(kind=dp),intent(in) :: theta
               integer, intent(in) :: l, m, n
               cosb=dcos(theta/2.d0)
               sinb=dsin(theta/2.d0)
               som=0.d0
               i=max(0,m-n)
               do while ((l+m+i).GE.0.and.(l-n-i).GE.0.and.(i+n+m).GE.0)
               som=som+APARITY(i)*(dsqrt(fac10(l+m)*fac10(l-m)*fac10(l+n)*fac10(l-n))/&
                       (fac10(l+m-i)*fac10(l-n-i)*fac10(i)*fac10(i+n-m)))&
                       *cosb**(2*l+m-n-2*i)*sinb**(2*i+n-m)
               i=i+1
               enddo
               dlmn=som
       return
       end function dlmn
 
 function fac10 (n)
  implicit none
! -----------------------------------------------
! function fac10(n) calculates factorial(n)/10**n
! -----------------------------------------------
! input: integer n >= 0 (you may want to check this
!        in the program calling this function)
! -----------------------------------------------
  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!------------------------------------------------
!      formal arguments
!------------------------------------------------
  integer, intent(in) :: n
!------------------------------------------------
!      local variables
!------------------------------------------------
  integer :: i
  real(wp) :: fac10, q
! -----------------------------------------------
  if (n == 0) then
     fac10 = 1.0_wp
  else
     fac10 = 1.0_wp
     q = 1.0_wp
     do i = 1, n
        fac10 = fac10 * q / 10.0_wp
        q = q + 1.0_wp
     end do
  endif
!
  return
  end function fac10

