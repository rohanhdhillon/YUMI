subroutine imprime2(iprint,imp,n,K,text)
  implicit real(kind=8) (a-h,o-z)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: dp = 8
  real(kind=dp) , intent(in), dimension(n,n)   :: K
  integer  ,  intent(in)     :: iprint, imp, n
  character(len=10), intent(in) :: text
  integer :: irest,i,j
  write(imp,*) text
  irest=(n/10)
  do j=1,irest
     write(imp,'(10I14)') (i,i=(j-1)*10+1,10*j)
     do i=1,n
        write(imp,'(I3,10E14.6)') i,K(i,(j-1)*10+1:j*10)
     enddo
  enddo
  if (irest*10 .NE. n) then
     write(imp,'(10I14)') (i,i=(j-1)*10+1,n)
     do i=1,n
        write(imp,'(I3,10E14.6)') i,K(i,irest*10+1:n)
     enddo
  endif
  return
end subroutine imprime2
