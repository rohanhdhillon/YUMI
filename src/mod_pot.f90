module mod_pot
        use mod_constantes, only :dp
  implicit real(kind=8) (a-h,o-z)
!  implicit none
!integer, parameter :: dp = 8
  real(kind=dp),allocatable, dimension(:) :: ri,hh,vvl
  real(kind=dp),allocatable, dimension(:,:) :: fi,fs
  real(kind=dp) ::  r
  real(kind=dp), allocatable :: flsp(:)
  integer, allocatable :: indexl(:), indexj(:)
contains



end module mod_pot
