module mod_constantes
implicit real(kind=8) (a-h,o-z)
!  implicit none
integer, parameter :: dp = 8
integer, parameter :: qp = 8
character(len=4) :: prognb
character(len=40)     :: nomfich
real(kind=dp), dimension(0:1) :: fpar=[1.d0,-1.d0]
integer, dimension(0:5) :: p6tab=[0,0,1,1,0,1]
real(kind=dp), parameter :: pi=3.141592653589793d0,&
                            FMPRT=1836.104D0,&                  ! unité de masse atomique
                            A02A2=2.8002852d-01,&               ! Angsetrum**2 <--- ua**2
                            BFCT=16.857630d0,&                  !
                            uacmm1=219474.63067d0               ! ua    ------->  cm-1
!FMPRT=1836.15152D0
!FMPRT=1822.8874D0
!                            FMPRT=1822.8874D0,&                 ! unité de masse atomique

contains

end module mod_constantes
