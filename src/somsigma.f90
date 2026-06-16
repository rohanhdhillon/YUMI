program somsigma
  use mod_constantes
  implicit real(kind=8) (a-h,o-z)
  real(kind=dp), allocatable, dimension(:,:) :: sigma1,sigma2,sigma3,sigma4
  character(len=24) :: nomfich
  character(len=8) ::cnopen
  integer :: nopen


  !**********************************************************************************
  !
  !..........................................  lecture du nom des fichiers de données
  !
  !**********************************************************************************
  call get_command_argument(1,nomfich)
  call get_command_argument(2,cnopen)
  !
  !write(*,*)nom
  imp=10
  iprint=1
  read(cnopen,'(I8)')nopen
  write(*,*) nopen,cnopen
  allocate(sigma1(nopen,nopen))
  allocate(sigma2(nopen,nopen))
  allocate(sigma3(nopen,nopen))
  allocate(sigma4(nopen,nopen))
  open (11,file=trim(nomfich)//"-1.sig",status="unknown",form="unformatted")
  open (12,file=trim(nomfich)//"-2.sig",status="unknown",form="unformatted")
  open (13,file=trim(nomfich)//"-3.sig",status="unknown",form="unformatted")
  open (14,file=trim(nomfich)//"-4.sig",status="unknown",form="unformatted")
  open (10,file=trim(nomfich)//".sig",status="unknown")
  read(11)sigma1
  read(12)sigma2
  read(13)sigma3
  read(14)sigma4
  sigma1=sigma1+sigma2+sigma3+sigma4
  call imprime2(iprint,imp,nopen,sigma1,'sigma     ')

end program somsigma
