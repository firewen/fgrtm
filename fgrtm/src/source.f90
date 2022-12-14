      subroutine SourceVector_D(kn,o)
!c     .......................

      include     "green_com.inc"
      complex*16  wru_0(2,2), wrd_0(2,2), wru_1(2,2), wrd_1(2,2)
      complex*16  wru_2(2,2), wrd_2(2,2)
      complex*16  csns, cpns, xc, vss, vps
      real*8      mus

!c     BASIC QUANTITIES:
!c     ................
!?      exd(1) = cdexp( -cpn(ls)*( z(ls)-zs ) )
!?      exd(2) = cdexp( -csn(ls)*( z(ls)-zs ) )
!?      exu(1) = cdexp( -cpn(ls)*(zs-z(ls-1)) )
!?      exu(2) = cdexp( -csn(ls)*(zs-z(ls-1)) )
	exu(1) = cdexp( -cpn(ls)*( z(ls)-zs ) )
      exu(2) = cdexp( -csn(ls)*( z(ls)-zs ) )
      exd(1) = cdexp( -cpn(ls)*(zs-z(ls-1)) )
      exd(2) = cdexp( -csn(ls)*(zs-z(ls-1)) )
!c     .....
      csns = csn(ls)
      cpns = cpn(ls)
      vss  = vs(ls)
      vps  = vp(ls)
      mus  = mu(ls)
      xc   = 2.*vps*cpns*csns*mus*o/vss

!c     SOURCE VECTORS:
!c      else if ( SourceType.eq.'D' ) then
!c        case 3: double-couple source:
!c        .......
      cigma0(1) = vss*cpns*cpns*csns/xc
      cigma0(2) = -kn*vps*cpns*csns/xc
      cigma1(1) =-2*vss*csns*cpns*kn/xc
      cigma1(2) = vps*(kn*kn+csns*csns)*cpns/xc
      cigma2(1) = vss*kn*kn*csns/xc
      cigma2(2) =-vps*cpns*csns*kn/xc
      cigma2_sh = kn/(2*mus*csns)
      cigma1_sh = 1.0/(2*mus)
      do j=1,2
        do i=1,2
         wru_0(i,j)=exu(i)*grdu(ls,i,j)+unit(i,j)
         wrd_0(i,j)=unit(i,j)+exd(i)*grud(ls-1,i,j)
         wru_1(i,j)=exu(i)*grdu(ls,i,j)-unit(i,j)
         wrd_1(i,j)=unit(i,j)-exd(i)*grud(ls-1,i,j)
        end do
      end do
      su0_2=(1d0+exu(2)*grdu0(ls))*cigma2_sh
      sd0_2=(exd(2)*grud0(ls-1)+1d0)*cigma2_sh
      su0_1=(1d0-exu(2)*grdu0(ls))*cigma1_sh
      sd0_1=(exd(2)*grud0(ls-1)-1d0)*cigma1_sh
      do i=1,2
         su_0(i) = wru_0(i,1)*cigma0(1) + wru_0(i,2)*cigma0(2)
         sd_0(i) = wrd_0(i,1)*cigma0(1) + wrd_0(i,2)*cigma0(2)
         su_1(i) = wru_1(i,1)*cigma1(1) + wru_1(i,2)*cigma1(2)
         sd_1(i) = wrd_1(i,1)*cigma1(1) + wrd_1(i,2)*cigma1(2)
         su_2(i) = wru_0(i,1)*cigma2(1) + wru_0(i,2)*cigma2(2)
         sd_2(i) = wrd_0(i,1)*cigma2(1) + wrd_0(i,2)*cigma2(2)
      end do         

      return
      end 
