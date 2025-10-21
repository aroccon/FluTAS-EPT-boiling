module mod_random
use mod_types
implicit none
private
public  random_value
contains
!
subroutine random_value(alow,ahigh,idum,iy,iv,outcome)
  implicit none
  !
  integer, intent(inout)       :: idum
  integer, intent(inout)        :: iy
  real(rp), intent(in) :: alow, ahigh
  real, dimension(1:32), intent(inout) :: iv
  real(rp)                , intent(out)  :: outcome
  integer :: j
  real :: iseed, ran1

  call randf_2(idum,iy,iv,ran1,iseed)
  idum = iseed
  outcome = (ran1*(ahigh-alow))+alow

  return
end subroutine random_value
!
subroutine randf_2(idum,iy,iv,ran1,iseed)
implicit none
integer, intent(inout) :: idum
integer, intent(inout) :: iy
real, dimension(1:32), intent(inout) :: iv
real, intent(out) :: iseed, ran1
integer, parameter :: ntab=32
real, parameter :: ia=16807., im=2147483647., am=1./im, iq=127773.
real, parameter :: ir=2836., ndiv=1.+(im-1.)/ntab, eps=1.2e-7, rnmx=1.-eps

integer :: j
real :: k

!! initial iseed (idum) is negative

if(idum .le. 0 .or. iy .eq. 0) then
   idum = max(-idum,1)
   do j=ntab+8,1,-1
     k=floor(idum/iq)
     idum=ia*(idum-k*iq)-ir*k
     if(idum .lt. 0.) then
       idum = idum + im
     endif
     if (j .le. ntab) then
       iv(j)=idum
     endif
   enddo
   iy = iv(1)
endif

k = floor(idum/iq)
idum = ia*(idum-k*iq)-ir*k

if (idum .le. 0) then
   idum = idum + im
endif

j=1+floor(iy/ndiv)
iy=iv(j)
iv(j)=idum
iseed=idum
ran1=min(am*iy,rnmx)

return
end subroutine randf_2
!
end module mod_random
