program yumi
  use mod_constantes
  use mod_fact
  use mod_input
  !$ use omp_lib
  implicit real(kind=dp) (a-h,o-z)
  integer :: input8,output,imp,ipot
  integer :: i, j, kk, l, ndim
  integer :: itype,nlev,lmin,lmax,jmax,lamax,jtot,nmax,nopen,nopen1
  integer :: itime(8),itime2(8)
  CHARACTER (LEN = 12) :: REAL_CLOCK (3)
  character(len=40)     :: label
  REAL time_begin, time_end
  namelist /input/ label,itype,rmax,rmin,npas,E,jtotl,jtotu,jtstep

  call factrl

1000 format('*              l1max =',I4,'        mmax =',I4,'        l2max =',I4,'         lmax =',I4,'               *')
1010 format('*                      Date      ',I4,'/',I2,'/',I4,'       Heure   ',I2,'h',I2,'mn                              *')
1020 format('*     Propateur de rmin =  ',F6.3,'   à rmax =  ',F6.3,'   avec un pas de =  ',I5,'                   * ')
1030 format('*     avec        jtotl =  ',I6,'  à jtotu =  ',I6,'   avec un pas de =  ',I5,'                   * ')
1005 format('*                        ',A40,'                  * ')
1040 format('*                                  Energie E = ',F10.4,'                                        * ')
2000 format(A98)
3000 format(A138)

  !      DATA ECONV/1.D0,0.6950387D0, 3.335640952D-5,3.335640952D-2,
  !     1   8065.5410D0,5.0341125D+15,219474.63067D0,
  !     2   83.593461D0,349.9891D0/

  !**********************************************************************************
  !
  !..........................................  lecture du nom des fichiers de données
  !
  !**********************************************************************************
  call get_command_argument(1,nomfich)
!  call get_command_argument(2,prognb)
  !
  write(*,*)trim(nomfich)

  input8=8
  output=9
  imp=output
  ipot=10
  itmat=15
  is11=11
  open (input8,file=trim(nomfich)//".inp",status="old")
  open (ipot,file=trim(nomfich)//".pot",status="old")
  open (output,file=trim(nomfich)//".out",status="unknown")

  label=' '
  itype=1
  rmin=2.d0
  rmax=30.d0
  npas=10
  E=100.d0
  jtotl=0
  jtotu=20
  jtstep=1
  read(input8,nml=input)
!  rmin=rmin/0.529177d0
!  rmax=rmax/0.529177d0
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                                         YUMI                                                   *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*                         (YUMI est l arc asymétrique japonais)                                 **'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*                           Programme  collision   moléculaire                                  **'
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                                                                                                *'
     write(output,1005) label
     write(output,2000)'*                                                                                                *'
     write(output,1020) rmin,rmax,npas
     write(output,2000)'*                                                                                                *'
     write(output,1030) jtotl,jtotu,jtstep
     write(output,2000)'*                                                                                                *'
     write(output,1040) E
     write(output,2000)'*                                                                                                *'

  time_begin = OMP_GET_WTIME()
  call date_and_time(real_clock(1),real_clock(2),real_clock(3),itime)

  select case(itype)
  case(1)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                             Atom − linear rigid rotor scattering                               *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*    A.M. Arthurs and A. Dalgarno, Proc. Roy. Soc. A256 540 (1963); doi:10.1098/rspa.1960.0125   *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'

     call flush(output)
!     call c0001(nomfich)
  case(3)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                          Rotateurs Linéaire -  Linéaire                                        *'
     write(output,2000)'*                                                                                                *'
!     write(output,2000)'*          Phillips et al. J. Chem. Phys. 101, 5824 (1994); doi: 10.1063/1.467297                *'
     write(output,2000)'*          Green    J. Chem. Phys. 62, 2271 (1975); doi: 10.1063/1.430752                        *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'

     call c0003(nomfich)
  case(4)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie Asymétrique -  Rotateur Linéaire                                 *'
     write(output,2000)'*                                                                                                *'
!     write(output,2000)'*          Heil     et al. J. Chem. Phys. 68, 2562 (1978); doi: 10.1063/1.436115                 *'
     write(output,2000)'*          Phillips et al. J. Chem. Phys. 102, 6024 (1995); doi: 10.1063/1.469337                *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'

     flush(output)
!     call c0004(nomfich)
  case(5)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie symétrique -   Atome sans structure                             *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*         Green   J. Chem. Phys. 64, 3463 (1976);    doi: 10.1063/1.432640                       *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'
     flush(output)

!     call c0005(nomfich)
  case(10)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie symétrique -  Rotateur Linéaire                                  *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*          Rist     et al. J. Chem. Phys. 98, 4662 (1993); doi: 10.1063/1.464970                 *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'
     flush(output)

!     call c0010(nomfich)
  case(11)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie symétrique -  Rotateur Linéaire                                  *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*                         Qianli Ma  27 august 2013  (CH3 - H2)                                  *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'
     flush(output)

!     call c0011(nomfich)
  case(12)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie symétrique -  Rotateur Linéaire                                  *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*               Alison R. Offer  Thèse Septembre 1990 Durham Université                          *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'
     flush(output)

!     call c0012(nomfich)
  case(15)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie symétrique -   Atome sans structure                             *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*         Green   J. Chem. Phys. 64, 3463 (1976);    doi: 10.1063/1.432640                       *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'
     flush(output)

!     call c0015(nomfich)
  case(40)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie Asymétrique -  Rotateur Linéaire                                 *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*          Phillips et al. J. Chem. Phys. 102, 6024 (1995); doi: 10.1063/1.469337                *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'

     flush(output)
!     call c0040(nomfich)
  case(140)
     write(output,2000)'**************************************************************************************************'
     write(output,2000)'*                        Toupie Asymétrique -  Rotateur Linéaire                                 *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'*          Phillips et al. J. Chem. Phys. 102, 6024 (1995); doi: 10.1063/1.469337                *'
     write(output,2000)'*                                                                                                *'
     write(output,2000)'**************************************************************************************************'
     write(output,1010)itime(3), itime(2), itime(1), itime(5), itime(6)
     write(output,2000)'**************************************************************************************************'

     flush(output)
!     call c0140(nomfich)
  case default
  end select
  write(output,2000)'**************************************************************************************************'
  time_end = OMP_GET_WTIME()
  call date_and_time(real_clock(1),real_clock(2),real_clock(3),itime2)
     write(output,1010)itime2(3), itime2(2), itime2(1), itime2(5), itime2(6)
     t1=itime2(1)-itime(1)   ! années
     t2=itime2(2)-itime(2)   ! mois
     t3=itime2(3)-itime(3)   ! jours
     t4=itime2(4)-itime(4)   !      
     t5=itime2(5)-itime(5)   ! heures
     t6=itime2(6)-itime(6)   ! minutes
     t7=itime2(7)-itime(7)   ! secondes
     t8=itime2(8)-itime(8)   ! 1/1000 de seconde
!  WRITE (output,*) t1,t2,t3,t4,t5,t6,t7,t8 
     t8=((t8/1000.)+t7)+(((t4*24.)+t5)*60.)*60.+t6*60.
  WRITE (output,'("         Wall time     ",F19.5,"      secondes")') t8
     write(output,2000)'**************************************************************************************************'
  WRITE (output,*) '        Temps CPU             ', time_end - time_begin, ' secondes'
  write(output,2000)'**************************************************************************************************'

end program yumi
