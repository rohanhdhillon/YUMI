subroutine pot
!implicit real(kind=8) (a-h,o-z)
use mod_constantes
use mod_pot
use mod_base
!$ use OMP_LIB
implicit real(kind=dp) (a-h,o-z)
!integer, parameter :: dp = 8
real(kind=dp), parameter   :: eps=0.0000001d0
real(kind=dp)   :: Ai,Bi,Ao,Bo,z
integer   :: n
integer        :: i,k,il
n=nmax
k=0


  if((r>=ri(0)-eps).and.(r<=ri(n)+eps)) then 
     do i=0,n-1
        if (r>=ri(i)) k=i
     enddo
  else
     !  print *,   'pot  :  revoir x !!!!',i,x(i),k
     if(r>=ri(n)+eps) then 
          !$OMP PARALLEL DO PRIVATE (il, Bo,Ao,z)
          do il=0,lamax
             Bo=dlog(dabs(fi(n-1,il)/fi(n,il)))/dlog(dabs(ri(n)/ri(n-1)))
             Ao=fi(n,il)*ri(n)**Bo
             z=Ao/(r**Bo)
             if(fi(n-1,il)*fi(n,il)<=0.d0.or.fi(n-1,il)/fi(n,il)<=1.0d0.or.Bo>8.0d1.or.Bo<=2.d0)then
                z=0.d0
                !write(2,*)r,il,z,Ao,Bo
             endif        
             vvl(il)=z
             !  z=z/100.d0
             !  z=0.d0
             !write(2,*)xp,z,Ao,Boo
          enddo
          !$OMP END PARALLEL DO
    else
    
      !$OMP PARALLEL DO PRIVATE (il,Bi,Ai,z)
      do il=0,lamax
         Bi=dlog(dabs(fi(0,il)/fi(1,il))/(ri(1)-ri(0)))
         Ai=fi(0,il)*dexp(Bi*ri(0))
         z=Ai*dexp(-Bi*r)
         if(fi(1,il)*fi(0,il)<=0.d0)z=0.d0
         vvl(il)=z
         !   z=f(0)
         !write(2,*)z,Ai,Bi
      enddo
      !$OMP END PARALLEL DO

    endif
    !  z=0.d0
    return
    !  stop
  endif

  if((r>=ri(0)-eps).and.(r<=ri(n)+eps)) then 
     !$OMP PARALLEL DO PRIVATE (il,z)
     do il=0,lamax
       z=fs(k,il)*((((ri(k+1)-r)**3)/(6.d0*hh(k)))-((hh(k)/6.d0)*(ri(k+1)-r)))
       z=z+fs(k+1,il)*((((r-ri(k))**3)/(6.d0*hh(k)))-((hh(k)/6.d0)*(r-ri(k))))
       z=z+(fi(k,il)/hh(k))*(ri(k+1)-r)
       z=z+(fi(k+1,il)/hh(k))*(r-ri(k))
       vvl(il)=z
   enddo
   !$OMP END PARALLEL DO
endif

!print *,   'l intervalle est:',k
!print *,   'pi(x)=',z
end subroutine pot
