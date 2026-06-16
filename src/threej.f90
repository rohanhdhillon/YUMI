       function threej(a,b,c,am,bm,cm)
!$ use OMP_LIB
               use mod_fact
!
!     program to compute the 3j symbol (a b c / am bm cm)
!     authors:  t. orlikowski and b. follmeg
!     current revision date: 4-may-1997
!
!   Adapted from xf3j routine in hibridon to replace thrj in molscat.
!   Faster and generally more accurate for large momenta than the
!   original thrj in molscat.
!
!   Suppressed all explicit reference to double precision intrinsics 
!   to support real*16 evaluation.  
!   Replaced real arguments by integer ones and changed to use integer
!   arithmetics whenever possible. Of course half-integer quantum numbers
!   are no more supported in this version. 
!
!   REAL*8 VERSION,   PV, 6 jan 2008.
! --------------------------------------------------------------------------
      implicit none
      real*8 threej
      real*8 x, delta
      real(kind=8) :: zero=0.d0, one=1.d0, half=0.5d0
      integer a,b,c,am,bm,cm, iabc, &
           iacbm,ibcam,iabmc,iamam,ibpbm,iapam,ibmbm,icpcm,icmcm, &
           minchi,maxchi, j1,j2,j3,j4, ll, a1,a2,a3,b1,b2,b3, &
           iaam, ibbm, xa, xb, ichi, l
!
      x = zero
! check for triangular conditions
      if ((c .gt. (a + b)) .or. (c .lt. abs(a - b))) goto 3
      iabc = a + b + c
      if ((am + bm + cm) .ne. 0) goto 3
      if (abs(am) .gt. a) goto 3
      if (abs(bm) .gt. b) goto 3
      if (abs(cm) .gt. c) goto 3
      if (am.eq.0 .and. bm.eq.0 .and.(mod(iabc,2)).ne.0) goto 3
      iacbm = a - c + bm
      ibcam = b - c - am
      iabmc = a + b - c
      iamam = a - am
      ibpbm = b + bm
      iapam = a + am + 1
      ibmbm = b - bm + 1
      icpcm = c + cm + 1
      icmcm = c - cm + 1
      minchi = max(0,ibcam,iacbm)
      maxchi = min(iabmc,iamam,ibpbm) - minchi
      iabmc = iabmc - minchi
      iaam = iamam - minchi
      ibbm = ibpbm - minchi
      ibcam = minchi - ibcam
      iacbm = minchi - iacbm
      iamam = iamam + 1
      ibpbm = ibpbm + 1
! compute delta
      j1 = iabc - (a+a) + 1
      j2 = iabc - (b+b) + 1
      j3 = iabc - (c+c) + 1
      j4 = iabc + 2
      delta = fact(j1-1) + fact(j2-1) + fact(j3-1) - fact(j4-1)
      ll=0
      x = one
      if (maxchi.le.0) goto 2
      a1 = ibcam + maxchi + 1
      a2 = iacbm + maxchi + 1
      a3 = minchi + maxchi + 1
      b1 = iabmc - maxchi
      b2 = iaam - maxchi
      b3 = ibbm - maxchi
      do ichi = 1, maxchi
         xa = (b1 + ichi) * (b2 + ichi) * (b3 + ichi)
         xb = ( (a1 - ichi) * (a2 - ichi) * (a3 - ichi) )
         x = one - (xa * x) / xb
      enddo
      if (x.eq.zero) then
         goto 3
      else if (x.lt.zero) then
         x = -x
         ll = 1
      endif
 2    x = log(x) - fact(iabmc+1-1) - fact(iaam+1-1) - fact(ibbm+1-1) - fact(ibcam+1-1) &
          - fact(iacbm+1-1) - fact(minchi+1-1)
      x = x+x + fact(iapam-1) + fact(iamam-1) + fact(ibpbm-1) + fact(ibmbm-1) + &
           fact(icpcm-1) + fact(icmcm-1) + delta
      x = exp( x * half )
      l = ll + minchi + b - a + cm
      if( mod(l,2) .ne. 0) x = -x
 3    threej = x
      return
      end function threej
