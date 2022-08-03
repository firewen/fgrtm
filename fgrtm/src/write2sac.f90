	subroutine write2sac(fname,dt,nt,yarray,tp,ts)
	implicit none
	include "sacf.h"

	integer :: SACMAX
	parameter (SACMAX=10000)

	character*20 fname
	integer :: nt
	real*8 :: dt,yarray(nt)
	real*8 :: tp,ts
	real :: delta,tmp(nt),t(nt)
	integer :: nerr
	real :: beg
	
	integer :: i

	beg = 0.0

	delta = sngl(dt)
	tmp = sngl(yarray)

	do i=1,nt
		t=(i-1)*delta
	end do
!	write(*,*)fname
	call newhdr()
	call setnhv('npts',nt,nerr)
	call setfhv('b',beg,nerr)
	call setfhv('e',t(nt),nerr)
	call setfhv('delta',delta,nerr)
	call setfhv('T0',sngl(tp),nerr)
	call setfhv('T1',sngl(ts),nerr)
	call wsac0(fname,t,tmp,nerr)
!	call wsac1(fname,tmp,nt,beg,delta,nerr)

	if (nerr.ne.0) then
		write(*,*)'Error writing SAC File:',fname
		call exit(-1)
	end if
	
	return
	end subroutine write2sac
