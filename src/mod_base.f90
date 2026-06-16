module mod_base
        use mod_constantes, only : dp
implicit real(kind=8) (a-h,o-z)
!  implicit none
!integer, parameter :: dp = 8
  integer :: lamax, nmax
  integer , allocatable, dimension(:) :: jtab, ltab,  l1,m1, l2, m2, l
  integer , allocatable, dimension(:) :: j1tab, katab, j2tab, j12tab, jjtab,ptab,ka2tab
  real(kind=dp) , allocatable, dimension(:) :: k2tab,ktab
  real(kind=dp) , allocatable, dimension(:,:) :: valp,epsk,alpha
  real(kind=dp) , allocatable, dimension(:,:) :: valp2,epsk2,alpha2
  real(kind=dp) , allocatable, dimension(:,:,:) :: vec
  real(kind=dp) , allocatable, dimension(:,:,:) :: vec2


contains



end module mod_base
