!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
      function w3js(j1,j2,j3,m1,m2,m3)
  use mod_fact
      integer :: m1, m2, m3, j1, j2, j3
            integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
!            integer :: iii(0:3)
            real(kind=qp) :: denom, cc, cc1, cc2

  real(kind=qp) :: w3js, w6js, w9js, som1
!            iii(0)=1
!            iii(1)=-1
!            iii(2)=-1
!            iii(3)=-1
            w3js = 0.d0
            if (m1+m2+m3 /= 0) goto 1000
            ia = j1 + j2
            if (j3 > ia) goto 1000
            ib = j1 - j2
            if (j3 < abs(ib)) goto 1000
            jsum = j3 + ia
            ic = j1 - m1
            id = j2 - m2

            if (abs(m1) > j1) goto 1000
                  if (abs(m2) > j2) goto 1000
                  if (abs(m3) > j3) goto 1000
                  ie = j3 - j2 + m1
                  im = j3 - j1 - m2
                  zmin = max0(0,-ie,-im)
                  ig = ia - j3
                  ih = j2 + m2
                  zmax = min0(ig,ih,ic)
                  cc1 = fact(ig/2)+fact((j3+ib)/2)+fact((j3-ib)/2)-fact((jsum+2)/2)
            cc2 = fact((j1+m1)/2)+fact(ic/2)+fact(ih/2)+fact(id/2)+fact((j3-m3)/2)+fact((j3+m3)/2)
                  cc = 0.d0
                  do z = zmin, zmax, 2
                  phase=1.d0
                  if (mod(z,4) /= 0) phase= -phase
!                  phase=iii(z-((z/4)*4))*phase
                        denom = fact(z/2)+fact((ig-z)/2)+fact((ic-z)/2)+fact((ih-z)/2)+&
                              fact((ie+z)/2)+fact((im+z)/2)
                        cc = cc + phase * dexp(-denom + (cc1+cc2)/2.d0)
                  enddo
            cc = cc !* dexp((cc1+cc2)/2.d0)
                  if (mod(ib-m3,4) /= 0) cc = -cc
                  w3js = cc
                  if (abs(w3js) < 1.d-8) w3js = 0.d0
                  if (abs(w3js) > 1.d30) write(*,*) w3js
1000             return
      end function w3js
            
!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
      function w6js(j1,j2,j3,l1,l2,l3)
            use mod_fact
            integer :: j1,j2,j3,l1,l2,l3
            integer :: ia, ib, ic, id, ie, iif, ig, ih, sum1, sum2, sum3, sum4
            integer :: w, wmin, wmax, ii, ij, ik            
            !integer :: iii(0:3)
            real(kind=qp) :: omega, denom, theta1, theta2, theta3, theta4, theta
            real(kind=qp) :: w3js, w6js, w9js, som1

            !iii(0)=1
            !iii(1)=-1
            !iii(2)=-1
            !iii(3)=-1
            w6js = 0.d0
            ia = j1 + j2
            if (ia < j3) goto 1000
              ib = j1 - j2
              if (abs(ib) > j3) goto 1000
                 ic = j1 + l2
                 if (ic < l3) goto 1000
                    id = j1 - l2
                    if (abs(id) > l3) goto 1000
                       ie = l1 + j2
                       if (ie < l3) goto 1000
                         iif = l1 - j2
                         if (abs(iif) > l3) goto 1000
                           ig = l1 + l2
                           if (ig < j3) goto 1000
                             ih = l1 - l2
                             if (abs(ih) > j3) goto 1000
                               sum1=ia + j3
                               sum2=ic + l3
                               sum3=ie + l3
                               sum4=ig + j3
                               wmin = max0(sum1, sum2, sum3, sum4)
                               ii = ia + ig
                               ij = j2 + j3 + l2 + l3
                               ik = j3 + j1 + l3 + l1
                               wmax = min0(ii,ij,ik)
                               theta1 = fact((ia-j3)/2)+fact((j3+ib)/2)+fact((j3-ib)/2)&
                               -fact(sum1/2+1)
                               theta2 = fact((ic-l3)/2)+fact((l3+id)/2)+fact((l3-id)/2)&
                               -fact(sum2/2+1)
                               theta3 = fact((ie-l3)/2)+fact((l3+iif)/2)+fact((l3-iiF)/2)&
                               -fact(sum3/2+1)
                               theta4 = fact((ig-j3)/2)+fact((j3+ih)/2)+fact((j3-ih)/2)&
                               -fact(sum4/2+1)
                               theta = theta1 + theta2 + theta3 + theta4
                               omega = 0.d0
                               do w = wmin, wmax, 2
                                 phase=1.d0
                                 if (mod(w,4) /= 0) phase= -phase
!                                phase=iii(mod(w,4)-((mod(w,4)/4)*4))*phase
                                       denom = fact((w-sum1)/2)+fact((w-sum2)/2)+fact((w-sum3)/2)&
                                               +fact((w-sum4)/2)+fact((ii-w)/2)+fact((ij-w)/2)&
                                               +fact((ik-w)/2)
                                        omega = omega + phase * dexp(fact(w/2+1) - denom + theta/2.d0)
                               enddo        
                               w6js = omega !* dexp(theta/2.d0)
                               if (abs(w6js) < 1.d-8) w6js = 0.d0
                               if (abs(w6js) > 1.d30) write(*,*) w6js
1000             return
      end function w6js

      function w9js(j1,j2,j3,j4,j5,j6,j7,j8,j9)
  use mod_fact
            integer :: j1,j2,j3,j4,j5,j6,j7,j8,j9
            integer :: i, kmin, kmax, k
            real(kind=qp) :: x, s, x1, x2, x3
  real(kind=qp) :: w3js, w6js, w9js, som1
            
            kmin = abs(j1-j9)
                  kmax = j1 + j9
                  i = abs(j4-j8)
                  if (i > kmin) kmin = i
                  i = j4 + j8
                  if (i < kmax) kmax = i
                  i = abs(j2-j6)
                  if (i > kmin) kmin = i
                  i = j2 + j6
                  if (i < kmax) kmax = i
                  x = 0.d0
                  do k = kmin, kmax, 2
                        s = 1.d0
                        if (mod(k,2) /= 0) s = -1.d0
                        x1 = w6js(j1,j9,k,j8,j4,j7)
                  x2 = w6js(j2,j6,k,j4,j8,j5)
                  x3 = w6js(j1,j9,k,j6,j2,j3)
                  x = x + s*x1*x2*x3*dfloat(k+1)
                  enddo
            w9js = x
            return
      end function w9js

      subroutine factrl
              use mod_fact
            integer :: i
            fact(0) = 0.d0
            do i=1,nfactmax
               fact(I) = fact(I-1) + dlog(dble(I))
            enddo
      END subroutine factrl
            
