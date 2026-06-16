subroutine Percival(jtot,nlev,fl)
        use mod_constantes
        use mod_fact
        use mod_pot
        use mod_base, only : lamax,jtab,ltab
  implicit real(kind=dp) (a-h,o-z)
  integer :: i, j, k, l, ndim, lam, irest
  integer, intent(in) :: jtot,nlev
  real(kind=dp),intent(out),dimension(nlev,nlev,0:lamax)::fl
  real(kind=dp) :: fac,som, f3j0
  real(kind=dp),allocatable, dimension(:) :: f
  real(kind=qp) :: w3js, w6js, w9js, som1
  eps=1.d-10
  fl(:,:,:)=0.d0
  do i=1,nlev
          do j=1,i
        fac=dsqrt(dfloat(2*jtab(i)+1)*dfloat(2*jtab(j)+1)*dfloat(2*ltab(i)+1)*dfloat(2*ltab(j)+1))
        fac=fac*(-1.d0)**(jtab(i)+jtab(j)-jtot)
        l1=max0(iabs(jtab(i)-jtab(j)),iabs(ltab(j)-ltab(i)))
        l2=min0(jtab(i)+jtab(j),ltab(j)+ltab(i))
        l3=min0(lamax,l2)
           do lam=l1,l3
              som1=threej0(jtab(j),jtab(i),lam)
              if(dabs(som1)<eps) cycle
              som2=threej0(ltab(j),ltab(i),lam)
              if(dabs(som2)<eps) cycle
              som3=sixj(jtab(i),ltab(i),jtot,ltab(j),jtab(j),lam)
              if(dabs(som3)<eps) cycle
              fl(i,j,lam)=som1*som2*som3*fac 
           enddo
     enddo
  enddo
    do lam=0,lamax
       do i=1,nlev
          do j=1,i
             fl(j,i,lam)=fl(i,j,lam)
          enddo
        enddo
     enddo
end subroutine Percival
