      function sixj(a,b,e,d,c,f)
!$ use OMP_LIB
              use mod_fact
!
!                                    | a  b  e |
!   program to compute the 6j symbol {         }
!                                    | d  c  f |
!   author: b. follmeg
!   current revision date: 4-may-1997
!
!   Modified version to replace sixj routine in Molscat. 
!   Note the transposition of c and e indices in the argument list
!   with respect to the 6j symbol. 
!   Suppressed all explicit reference to double precision intrinsics 
!   to support real*16 evaluation.  
!   Replaced real arguments by integer ones and changed to use integer
!   arithmetics whenever possible. Of course half-integer quantum numbers
!   are no more supported in this version. 
!
!   REAL*8 VERSION,   PV, 13 dec 2007.
! -------------------------------------------------------------------
      implicit none
      real*8 sixj
      real*8 x, delta
      real(kind=8) :: zero=0.d0, one=1.d0, half=0.5d0
      integer a,b,e,d,c,f
      integer iabe, idce, iacf, idbf
      integer iabdc, iaedf, ibecf, minchi, maxchi
      integer i2a, i2b, i2c, i2d, i2e, i2f, j1,j2,j3,j4
      integer abdc, aedf, becf, abe, dce, acf, dbf, xa, xb 
      integer ipower, ii, ichi
      x=zero
! check triangular conditions for triad ( a b e)
      if ((e .gt. (a + b)) .or. (e .lt. abs(a - b))) goto 40
      iabe = a + b + e
! check triangular conditions for triad ( d c e)
      if ((e .gt. (c + d)) .or. (e .lt. abs(c - d))) goto 40
      idce = d + c + e
! check triangular conditions for triad ( a c f)
      if ((f .gt. (a + c)) .or. (f .lt. abs(a - c))) goto 40
      iacf = a + c + f
! check triangular conditions for triad ( d b f)
      if ((f .gt. (d + b)) .or. (f .lt. abs(d - b))) goto 40
      idbf = d + b + f
      iabdc = a + b + d + c
      iaedf = a + e + d + f
      ibecf = b + e + c + f
      minchi = max(iabe,idce,iacf,idbf,-1)
      maxchi = min(iabdc,iaedf,ibecf) - minchi
! indices for deltas
      delta = zero
      i2a = a+a - 1
      i2b = b+b - 1
      i2c = c+c - 1
      i2d = d+d - 1
      i2e = e+e - 1
      i2f = f+f - 1
! delta(abe)
      j1 = iabe - i2a
      j2 = iabe - i2b
      j3 = iabe - i2e
      j4 = iabe + 2
      delta = delta + fact(j1-1) + fact(j2-1) + fact(j3-1) - fact(j4-1)
! delta(dce)
      j1 = idce - i2d
      j2 = idce - i2c
      j3 = idce - i2e
      j4 = idce + 2
      delta = delta + fact(j1-1) + fact(j2-1) + fact(j3-1) - fact(j4-1)
! delta(acf)
      j1 = iacf - i2a
      j2 = iacf - i2c
      j3 = iacf - i2f
      j4 = iacf + 2
      delta = delta + fact(j1-1) + fact(j2-1) + fact(j3-1) - fact(j4-1)
! delta(dbf)
      j1 = idbf - i2d
      j2 = idbf - i2b
      j3 = idbf - i2f
      j4 = idbf + 2
      delta = delta + fact(j1-1) + fact(j2-1) + fact(j3-1) - fact(j4-1)
      delta = half * delta
      iabdc = iabdc - minchi
      iaedf = iaedf - minchi
      ibecf = ibecf - minchi
      iabe = minchi - iabe
      idce = minchi - idce
      iacf = minchi - iacf
      idbf = minchi - idbf
      abdc = iabdc - maxchi
      aedf = iaedf - maxchi
      becf = ibecf - maxchi
      abe = maxchi + iabe + 1
      dce = maxchi + idce + 1
      acf = maxchi + iacf + 1
      dbf = maxchi + idbf + 1
! loop over chi
      x = one
      ipower = 0
      if (maxchi .le. 0) goto 30
      ii = minchi + maxchi + 2
      do ichi = 1, maxchi
         xa = ( (abdc + ichi) * (aedf + ichi) * (becf + ichi) ) &
              * (ii - ichi)
         xb = (abe - ichi) * (dce - ichi) * (acf - ichi) * (dbf - ichi)
         x = one - ( xa * x ) / xb
      enddo
      if (x.eq.zero) then
         goto 40
      else if (x.lt.zero) then
         x = -x
         ipower = 1
      endif
 30   x = log(x) + fact(minchi+2-1)-fact(iabdc+1-1)-fact(iaedf+1-1)-fact(ibecf+1-1) &
                  - fact(iabe+1-1)-fact(idce+1-1)-fact(iacf+1-1)-fact(idbf+1-1) + delta
      ipower = ipower + minchi
      x = ((-1)**ipower) * exp(x)
 40   sixj = x
      return
      end function sixj
