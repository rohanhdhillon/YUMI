      function threej0(a,b,c)
              use mod_fact
!
!     program to compute the 3j symbol (a b c / 0 0 0)
!     author:  b. follmeg
!     current revision date: 14-dec-87
!     revised by mha 4-may-1997
!
!     adapted from xf3jm0 routine in hibridon to replace threej in molscat.
!     much faster than original threej in molscat (no loops). 
!
!     REAL*8 VERSION,   PV, 13 dec 2007.
! --------------------------------------------------------------------------
      implicit none
      real*8 threej0
      real*8  x, delta
      real(kind=8) :: zero=0.d0,half=0.5d0
      integer a, b, c, iabc, ig, j1,j2,j3,j4, ip
!
      x = zero
! check for triangular conditions
      if ((c .gt. (a + b)) .or. (c .lt. abs(a - b))) goto 100
      iabc = a + b + c
! check for even sum
      ig = iabc / 2
      if (ig*2 .ne. iabc) goto 100
! compute delta
      j1 = iabc - (a+a) + 1
      j2 = iabc - (b+b) + 1
      j3 = iabc - (c+c) + 1
      j4 = iabc + 2
      delta = half * (fact(j1-1) + fact(j2-1) + fact(j3-1) - fact(j4-1))
      j1 = ig + 1
      j2 = ig - a + 1
      j3 = ig - b + 1
      j4 = ig - c + 1
      x = delta + fact(j1-1) - fact(j2-1) - fact(j3-1) - fact(j4-1)
      x = exp(x)
      ip = ig + c + a - b
      if( mod(ip,2) .ne. 0) x = -x
 100  threej0 = x
      return
      end function threej0
