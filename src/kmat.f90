subroutine kmat(RP,nbv,ltab,ktab,rmax,K)
integer, intent(in) :: nbv
integer, intent(in), dimension(nbv) :: ltab
real(kind=8), intent(in), dimension(nbv) :: ktab
real(kind=8), intent(in), dimension(nbv,nbv) ::RP
real(kind=8), intent(in) :: rmax
real(kind=8), intent(out), dimension(nbv,nbv) :: K
real(kind=8), allocatable,  dimension(:,:) :: Jn,Nn,Jnp,Nnp,gg
real(kind=8), allocatable, dimension(:) :: work
real(kind=8) :: z
integer ::info,lwork,irest
integer, allocatable, dimension(:) :: ipiv
allocate(ipiv(nbv))
allocate(Jn(nbv,nbv),Nn(nbv,nbv),Jnp(nbv,nbv),Nnp(nbv,nbv),gg(nbv,nbv))
lwork=64*nbv
allocate(work(lwork))
Jn=0.d0
Jnp=0.d0
Nn=0.d0
Nnp=0.d0
do i=1,nbv
if(ktab(i)>=0.d0)then
z=rmax*dsqrt(ktab(i))
call bessel(ltab(i),z,Jn(i,i),Jnp(i,i),Nn(i,i),Nnp(i,i))
else
Jn(i,i)=1.d0
Nn(i,i)=1.d0
Jnp(i,i)=1.d0
Nnp(i,i)=1.d0
endif
enddo
call dgemm('N','N',nbv,nbv,nbv,-1.d00,RP,nbv,Jn,nbv,1.d00,Jnp,nbv)
!  irest=(nbv/10)
!    do j=1,irest
!      print *, (j-1)*10+1
!      do i=1,nbv
!       write(*,'(I3,10E14.6)') i,Jnp(i,(j-1)*10+1:j*10)
!      enddo
!    enddo
!    if (irest*10 .NE. nbv) then
!      print *, (j-1)*10+1
!      do i=1,nbv
!        write(*,'(I3,10E14.6)') i,Jnp(i,irest*10+1:nbv)
!      enddo
!    endif
call dgemm('N','N',nbv,nbv,nbv,1.d00,RP,nbv,Nn,nbv,-1.d00,Nnp,nbv)
call dgetrf(nbv,nbv,Nnp,nbv,ipiv,info)
call dgetri(nbv,Nnp,nbv,ipiv,work,lwork,info)
call dgemm('N','N',nbv,nbv,nbv,1.d00,Jnp,nbv,Nnp,nbv,0.d00,K,nbv)
end subroutine kmat
