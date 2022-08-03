!c
!c	"GRT_COEFS.F" is used to calculate the generalized R/T coeffecients
!c
	subroutine grt_coefs(kn,o)
	
	include	"green_com.inc"
	
	complex*16	a(4,4),b(4,4),ab(4,4),eu(4,4),ed(4,4),cu,cd
	complex*16	tempa,tempb
	real*8		hu,hd

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc      Calculate generalized R/T coefficients(j=1,2,...,N)  ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!c	grud at the free-surface (lay=0)
!c
	if (ls.ne.1) then
		hu=z(1)-z(0)
		exu(1)=cdexp(-cpn(1)*hu)
		exu(2)=cdexp(-csn(1)*hu)
	else
		exu(1) = cdexp(-cpn(1)*(zs-z(0)))
		exu(2) = cdexp(-csn(1)*(zs-z(0)))
	end if

!c	SH waves
	grud0(0)=exu(2)

!c	P-SV waves
	call mtxe(1, kn, o, eu)
	do j=1,2
		do i=1,2
			a(i,j) = eu(i+2,j)
			b(i,j) = eu(i+2,j+2)
		end do
	end do
	
	call mtxinv_2(a)

	do j=1,2
		do i=1,2
			grud(0,i,j) = -(a(i,1)*b(1,j)+a(i,2)*b(2,j))*exu(j)
		end do
	end do

!c	grud an gtu for lay=1,2,...,(ls-2)-th interfaces:

	do lay=1,ls-1
		hd = z(lay)-z(lay-1)
		exd(1) = cdexp(-cpn(lay)*hd)
		exd(2) = cdexp(-csn(lay)*hd)
!		exd(1) = exu(1)
!		exd(2) = exu(2)
		if (lay.ne.ls-1) then
			hu = z(lay+1)-z(lay)
			exu(1) = cdexp(-cpn(lay+1)*hu)
			exu(2) = cdexp(-csn(lay+1)*hu)
		else
!c specifically deal with	grud an gtu for lay=(ls-1)-th interfaces:
			exu(1) = cdexp(-cpn(lay+1)*(zs-z(lay)))
			exu(2) = cdexp(-csn(lay+1)*(zs-z(lay)))
		end if
!		ex(1,1) = exu(1)
!		ex(1,2) = exu(2) 

!c	SH waves
		tempa = 0.5+mu(lay)*csn(lay)/(2d0*mu(lay+1)*csn(lay+1))
		tempb = 1d0-tempa

		gtu0(lay) = exu(2)/(tempa+tempb*exd(2)*grud0(lay-1))
		grud0(lay) = (tempb+tempa*exd(2)*grud0(lay-1))*gtu0(lay)

!c	P-SV waves
		call mtxe_i_4(lay+1, kn, o, a)
		call mtxe(lay, kn, o, b)
		call mtxmtp(4, a, b, ab)

		do j=1,2
			do i=1,2
				a(i,j) = ab(i+2,1)*exd(1)*grud(lay-1,1,j)	&
     						+ab(i+2,2)*exd(2)*grud(lay-1,2,j)	&
     						+ab(i+2,j+2)
				b(i,j) = ab(i,1)*exd(1)*grud(lay-1,1,j)		&
     						+ab(i,2)*exd(2)*grud(lay-1,2,j)	&
     						+ab(i,j+2)
			end do
		end do

		call mtxinv_2(a)

		do j=1,2
			do i=1,2
				gtu(lay,i,j) = a(i,j)*exu(j)
			end do
		end do
		do j=1,2
		   do i=1,2
			 grud(lay,i,j) = b(i,1)*gtu(lay,1,j)+b(i,2)*gtu(lay,2,j)
		   end do
		end do
	end do

!c END FOR grud an gtu for lay=1,2,...,(ls-1)-th interfaces:

!c	grdu an gtd for lay=(nly-1)-th interfaces:	
!c	because grdu(nly)=0
	lay = nly-1				!z(N) interfaces
	if (lay.ne.ls) then
		hd = z(lay)-z(lay-1)
		exd(1) = cdexp(-cpn(lay)*hd)
		exd(2) = cdexp(-csn(lay)*hd)
	else
		exd(1) = cdexp(-cpn(lay)*(z(ls)-zs))
		exd(2) = cdexp(-csn(lay)*(z(ls)-zs))
	end if

!c	SH waves	
	tempa = 0.5+mu(lay+1)*csn(lay+1)/(2d0*mu(lay)*csn(lay))
	tempb = 1d0-tempa

	gtd0(lay) = exd(2)/tempa
	grdu0(lay) = tempb*gtd0(lay)

!c	P-SV waves
	call mtxe_i_4(lay, kn, o, a)
	call mtxe(lay+1, kn, o, b)
	call mtxmtp(4, a, b, ab)

!	grud(N+1,i,j)=0
	do j=1,2
		do i=1,2
			a(i,j) = ab(i,j)
!     &				+ab(i,3)*exu(1)*grud(1,j)
!     &				+ab(i,4)*exu(2)*grud(2,j)
			b(i,j) = ab(i+2,j)
!     &				+ab(i+2,3)*exu(1)*grud(1,j)
!     &				+ab(i+2,4)*exu(2)*grud(2,j)
		end do
	end do

	call mtxinv_2(a)

	do j=1,2
		do i=1,2
			gtd(lay,i,j) = a(i,j)*exd(j)
		end do
	end do
	do j=1,2
	   do i=1,2
		 grdu(lay,i,j) = b(i,1)*gtd(lay,1,j)+b(i,2)*gtd(lay,2,j)
	   end do
	end do

!c	grdu an gtd for lay=lay=nly-2,...,ls-th interfaces:
	do lay=nly-2,ls,-1
		hu = z(lay+1)-z(lay)
		exu(1) = cdexp(-cpn(lay+1)*hu)
		exu(2) = cdexp(-csn(lay+1)*hu)
!		exu(1) = exd(1)
!		exu(2) = exd(2)
		if (lay.ne.ls) then
			hd = z(lay)-z(lay-1)
			exd(1) = cdexp(-cpn(lay)*hd)
			exd(2) = cdexp(-csn(lay)*hd)
		else
!c specifically deal with grdu an gtd for lay=ls-th interfaces:
			exd(1) = cdexp(-cpn(lay)*(z(lay)-zs))
			exd(2) = cdexp(-csn(lay)*(z(lay)-zs))
		end if

!c	SH waves
		tempa = 0.5+mu(lay+1)*csn(lay+1)/(2d0*mu(lay)*csn(lay))
		tempb = 1d0-tempa

		gtd0(lay) = exd(2)/(tempa+tempb*exu(2)*grdu0(lay+1))
		grdu0(lay) = (tempb+tempa*exu(2)*grdu0(lay+1))*gtd0(lay)

!c	P-SV waves
		call mtxe_i_4(lay, kn, o, a)
		call mtxe(lay+1, kn, o, b)
		call mtxmtp(4, a, b, ab)

		do j=1,2
			do i=1,2
				a(i,j) = ab(i,j)						&
     					+ab(i,3)*exu(1)*grdu(lay+1,1,j)		&
     					+ab(i,4)*exu(2)*grdu(lay+1,2,j)
				b(i,j) = ab(i+2,j)					&
     					+ab(i+2,3)*exu(1)*grdu(lay+1,1,j)         &
     					+ab(i+2,4)*exu(2)*grdu(lay+1,2,j)
			end do
		end do

		call mtxinv_2(a)

		do j=1,2
			do i=1,2
				gtd(lay,i,j) = a(i,j)*exd(j)
			end do
		end do
		do j=1,2
		   do i=1,2
			 grdu(lay,i,j) = b(i,1)*gtd(lay,1,j)+b(i,2)*gtd(lay,2,j)
		   end do
		end do
	end do

!c END FOR grdy an gtd for lay=nly-2,...,ls-th interfaces:

	return
	end