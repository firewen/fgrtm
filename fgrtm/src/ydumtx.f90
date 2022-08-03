
	subroutine ydumtx(e)

	include		"green_com.inc"
	complex*16	yu_psv(4,4),yd_psv(4,4),y_sh,temp1(4,4),temp2(4,4)
	complex*16	e(4,4),exdd(2),exuu(2),exds(2),exus(2)
	complex*16	u0,d0,u(2,2),d(2,2)

	exd(1) = cdexp(-cpn(ls)*(z(ls)-zs))
	exd(2) = cdexp(-csn(ls)*(z(ls)-zs))
	exu(1) = cdexp(-cpn(ls)*(zs-z(ls-1)))
	exu(2) = cdexp(-csn(ls)*(zs-z(ls-1)))
!?	exu(1) = cdexp(-cpn(ls)*(z(ls)-zs))
!?	exu(2) = cdexp(-csn(ls)*(z(ls)-zs))
!?	exd(1) = cdexp(-cpn(ls)*(zs-z(ls-1)))
!?	exd(2) = cdexp(-csn(ls)*(zs-z(ls-1)))

	exdd(1) = cdexp(-cpn(lo)*(z0-z(lo-1)))
      exdd(2) = cdexp(-csn(lo)*(z0-z(lo-1)))
      exuu(1) = cdexp(-cpn(lo)*(z(lo)-z0))
      exuu(2) = cdexp(-csn(lo)*(z(lo)-z0))

! prepare yu_psv(i,j), yd_psv(i,j) and y_sh
	do j=1,2
		do i=1,2
			temp1(i,j)=exd(i)*grdu(ls,i,j)
			temp2(i,j)=exu(i)*grud(ls-1,i,j)
		end do
	end do

	call mtxmtp(2,temp1,temp2,yu_psv)
	call mtxmtp(2,temp2,temp1,yd_psv)

	do j=1,2
		do i=1,2
			yu_psv(i,j)=unit(i,j)-yu_psv(i,j)
			yd_psv(i,j)=unit(i,j)-yd_psv(i,j)
		end do
	end do

      call mtxinv_2(yu_psv)
      call mtxinv_2(yd_psv)

	y_sh = 1d0/(1d0-exd(2)*grdu0(ls)*exu(2)*grud0(ls-1))

!cccccccccccccccc the case of j=s ccccccccccccccccccccccccccccc
	if (lo.eq.ls.and.z0.le.zs) then
! see (3-56d)
!?		exds(1) = cdexp(-cpn(ls)*(z0-zs))
!?		exds(2) = cdexp(-csn(ls)*(z0-zs))
! see (3-56c)
!		exus(1) = cdexp(-cpn(ls)*(zs-z0))
!		exus(1) = cdexp(-csn(ls)*(zs-z0))
		exds(1)=cdexp(-cpn(ls)*(z0-z(ls-1)))
		exds(2)=cdexp(-csn(ls)*(z0-z(ls-1)))
		exus(1)=cdexp(-cpn(ls)*(zs-z0))
		exus(2)=cdexp(-csn(ls)*(zs-z0))

! for SH waves
		u0=exds(2)*grud0(ls-1)+exus(2)
! for PSV waves
		do j=1,2
			do i=1,2
				u(i,j)=e(i,1)*exds(1)*grud(ls-1,1,j) &
     				+e(i,2)*exds(2)*grud(ls-1,2,j)		 &
     				+e(i,j+2)*exus(j)
			end do
		end do

		yu0=u0*y_sh
		do j=1,2
			do i=1,2
				yu(i,j)=u(i,1)*yu_psv(1,j)+u(i,2)*yu_psv(2,j)
			end do
		end do
!	write(*,*)"break"
	else if(lo.eq.ls.and.z0.gt.zs) then

		exds(1)=cdexp(-cpn(ls)*(z0-zs))
		exds(2)=cdexp(-csn(ls)*(z0-zs))
		exus(1)=cdexp(-cpn(ls)*(z(ls)-z0))
		exus(2)=cdexp(-csn(ls)*(z(ls)-z0))

		d0=exds(2)+exus(2)*grdu0(ls)
		do j=1,2
			do i=1,2
				d(i,j)=e(i,j)*exds(j)				&
     					+e(i,3)*exus(1)*grdu(ls,1,j)	&
     					+e(i,4)*exus(2)*grdu(ls,2,j)	
			end do
		end do

		yd0=d0*y_sh
		do j=1,2
			do i=1,2
				yd(i,j)=d(i,1)*yd_psv(1,j)+d(i,2)*yd_psv(2,j)
			end do
		end do

!ccccccccccccccc the case of j<s ccccccccccccccccccccc
	else if(lo.lt.ls) then
		u0=exdd(2)*grud0(lo-1)+exuu(2)
		do j=1,2
			do i=1,2
				u(i,j)=e(i,1)*exdd(1)*grud(lo-1,1,j)	&
     					+e(i,2)*exdd(2)*grud(lo-1,2,j)		&
    					+e(i,j+2)*exuu(j)
			end do
		end do

		do lay=ls-1,lo,-1
			if (lay.eq.ls-1) then
				yu0=gtu0(lay)*y_sh
				do j=1,2
					do i=1,2
						yu(i,j)=gtu(lay,i,1)*yu_psv(1,j)	&
     							+gtu(lay,i,2)*yu_psv(2,j)
					end do
				end do
			else
				yu0=gtu0(lay)*yu0
				do j=1,2
					do i=1,2
						temp1(i,j)=gtu(lay,i,1)*yu(1,j)		&
								+gtu(lay,i,2)*yu(2,j)
					end do
				end do
				do j=1,2
					do i=1,2
						yu(i,j)=temp1(i,j)
					end do
				end do
			end if
		end do

		yu0=u0*yu0
		do j=1,2
			do i=1,2
				temp1(i,j)=u(i,1)*yu(1,j)+u(i,2)*yu(2,j)
			end do
		end do
		do j=1,2
			do i=1,2
				yu(i,j)=temp1(i,j)
			end do
		end do

!ccccccccccccccccc the case of j>s cccccccccccccccccccccccc
	else if(lo.gt.ls) then
		d0=exdd(2)+exuu(2)*grdu0(lo)
		do j=1,2
			do i=1,2
				d(i,j)=e(i,j)*exdd(j)				&
     					+e(i,3)*exuu(1)*grdu(lo,1,j)	&
					+e(i,4)*exuu(2)*grdu(lo,2,j)
			end do
		end do

		do lay=ls+1,lo                               
			if (lay.eq.ls+1) then
				yd0=gtd0(ls)*y_sh
				do j=1,2
					do i=1,2           
						yd(i,j)=gtd(ls,i,1)*yd_psv(1,j)		&
     							+gtd(ls,i,2)*yd_psv(2,j)
					end do
				end do
			else
				yd0=gtd0(lay-1)*yd0
				do j=1,2
					do i=1,2
						temp1(i,j)=gtd(lay-1,i,1)*yd(1,j)	&
							+gtd(lay-1,i,2)*yd(2,j)
				end do
				end do
				do j=1,2
					do i=1,2
						yd(i,j)=temp1(i,j)
                      end do
				end do
			end if
		end do

		yd0=d0*yd0
		do j=1,2
			do i=1,2
				temp1(i,j)=d(i,1)*yd(1,j)+d(i,2)*yd(2,j)
			end do
		end do
		do j=1,2
			do i=1,2
				yd(i,j)=temp1(i,j)
			end do
		end do

	end if

	return
	end