      real*8 function ninej(a,b,c,d,e,f,g,h,i)
!$ use OMP_LIB
              use mod_fact
!
!                                        | a b c |
!... function to compute the 9j - symbol { d e f }
!                                        | g h i |
!
!    this is done by contracting three 6j - symbols.
!    (see edmonds, p. 101, eq. (6.4.3))
!     current revision date: 3-apr-91 by mha
!
!   Modified version to replace xninej routine in Molscat. 
!   Suppressed all explicit reference to double precision intrinsics 
!   to support real*16 evaluation.  
!   Replaced real arguments by integer ones and changed to use integer
!   arithmetics whenever possible. Of course half-integer quantum numbers
!   are no more supported in this version. 
!
!   REAL*8 VERSION,   PV, 13 dec 2007.
! -------------------------------------------------------------------
      implicit none
      real*8  x, delta, t1,t2,t3, sixj
      real(kind=8) :: zero=0.d0, one=1.d0, half=0.5d0
      integer a,b,c,d,e,f,g,h,i, kj1(6),kj2(6), kmx,kmn, sum,dif, is, kap, inorm
!
      x=zero
! check triangular conditions for triad ( a b c)
      if ((c .gt. (a + b)) .or. (c .lt. abs(a - b))) goto 150
! check triangular conditions for triad ( d e f)
      if ((f .gt. (d + e)) .or. (f .lt. abs(d - e))) goto 150
! check triangular conditions for triad ( g h i)
      if ((i .gt. (g + h)) .or. (i .lt. abs(g-h))) goto 150
! check triangular conditions for triad ( a d g)
      if ((g .gt. (a + d)) .or. (g .lt. abs(a - d))) goto 150
! check triangular conditions for triad ( b e h)
      if ((h .gt. (b + e)) .or. (h .lt. abs(b - e))) goto 150
! check triangular conditions for triad ( c f i)
      if ((i .gt. (c + f)) .or. (i .lt. abs(c - f))) goto 150
! restrict sum
! j1 and j2 are pairs that are coupled with kappa in triangular relations
      kj1(1)=a
      kj1(2)=h
      kj1(3)=b
      kj1(4)=d
      kj1(5)=a
      kj1(6)=f
      kj2(1)=i
      kj2(2)=d
      kj2(3)=f
      kj2(4)=h
      kj2(5)=i
      kj2(6)=b
      kmx=1000000000
      kmn=-1
      do is=1,6
         sum=kj1(is)+kj2(is)
         dif=abs(kj1(is)-kj2(is))
         if (sum .lt. kmx) kmx=sum
         if (dif .gt. kmn) kmn=dif
      enddo
      if (kmn .gt. kmx) then
        print *,' kmn = ',kmn,' .gt. kmx = ',kmx,' in xninej; abort'
        stop
      endif
! main sum
      do kap=kmn,kmx
         inorm=((-1)**(2*kap)) * (2*kap+1)
         t1=sixj(a,d,g,h,i,kap)
         t2=sixj(b,e,h,d,kap,f)
         t3=sixj(c,f,i,kap,a,b)
!         t1=sixj(a,d,i,h,g,kap)
!         t2=sixj(b,e,kap,d,h,f)
!         t3=sixj(c,f,a,kap,i,b)
         x=x+inorm*t1*t2*t3
      enddo
150   ninej=x
      return
      end function ninej
