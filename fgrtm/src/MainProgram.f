
      include   "green_com.inc"
      
      integer   fff1, fff2
      real*8    f_1, f_2, amp0 ,win ,freal, window, t, ifftcoef
      complex*16   st0, sw, rrr(n), fff(n), zzz(n)
      complex*16   integf(n/2,10), integt(n,10)
	real*8    tmp(n+5)
	real*8    urr(n), uff(n), uzz(n)
	character*20  fname

	real*8    tp, ts

!c  (0). Basic constants:
!c  -----------------------------
      aj = cmplx(0d0, 1d0)
      pi = 4d0*atan(1d0)
      pi2= 2d0*pi

!c  (1). Reading & Checking input:
!c  -------------------------------   
      call green_input

!c  (2). Basic parameters: 
!c  ----------------------
      call green_basic
!c
!c  (3). Radiation patterns (fai vs. str, dip & rak):
!c  -------------------------------------------------
	if (Mtype.eq.0) then
            f_1   = fai0 - str
            f_2   = 2d0*f_1
              fpsv01 = - sin(dip)*cos(dip)*sin(rake)
            fpsv02 = - 2d0*fpsv01
            fpsv1  = cos(2d0*dip)*sin(rake)*sin(f_1)
     &             - cos(  dip  )*cos(rake)*cos(f_1)
            fpsv2  = sin(dip)*( cos(rake)*sin(f_2)
     &             + cos(dip)*sin(rake)*cos(f_2) )
            fsh1   = - cos(2d0*dip)*sin(rake)*cos(f_1)
     &             - cos(  dip  )*cos(rake)*sin(f_1)
            fsh2   = sin(dip)*( cos(rake)*cos(f_2)
     &             - cos(dip)*sin(rake)*sin(f_2) )
      else if (Mtype.eq.1) then
              fpsv01 = (mxx+myy)/2
              fpsv02 = mzz
              fpsv1  = mxz*cos(fai0)+myz*sin(fai0)
              fpsv2  = -(myy-mxx)*cos(2*fai0)/2+mxy*sin(2*fai0)
              fsh1   = mxz*sin(fai0)-myz*cos(fai0)
              fsh2   = (myy-mxx)*sin(2*fai0)/2+mxy*cos(2*fai0)
      end if
!        f_1   = fai0 
!        f_2   = 2d0*f_1
!        fpsv01 = (mxx+myy)/2
!        fpsv02 = mzz
!        fpsv1  = mxz*cos(fai0)+myz*sin(fai0)
!        fpsv2  = -(myy-mxx)*cos(2*fai0)/2+mxy*sin(2*fai0)
!        fsh1   = mxz*sin(fai0)-myz*cos(fai0)
!        fsh2   = (myy-mxx)*sin(2*fai0)/2+mxy*cos(2*fai0)

      amp0 = M0

!c  (4). Calculate displacement spectra: Ur(r0,z0,fai0;o),
        nf2=m1+1
        nf1=1
         
        if (WinSwitch.eq.'OFF') then
            fff2 = m1+1 
            fff1 = 1   
        else 
            fff2 = min(int(f4/df), m1)+1
            fff1 = int(f1/df)+1
        end if 

        call green_spectra_r(fff1, fff2, integf)


!c      pre-FFT process:  
        do i = nf1,nf2
           freal = df*(i-1)         
           o     = pi2*freal - aj*oi
           st0   = sw(o)
           win   = window(freal)
	     ifftcoef = df*mt/pi2
           integt(i,:) = amp0*st0*win*integf(i,:)*ifftcoef
        end do

        do i=nf2+1, mt
           integt(i,:) = conjg(integt(mt+2-i,:))
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! directly output spectra
        rrr=integt(:,3)*fpsv2+integt(:,4)*fpsv1+integt(:,5)*fpsv01
     &      +integt(:,6)*fpsv02
        fff=integt(:,1)*fsh1+integt(:,2)*fsh2
        zzz=integt(:,7)*fpsv2+integt(:,8)*fpsv1+integt(:,9)*fpsv01
     &      +integt(:,10)*fpsv02
!        open(21,file=trim(Outname)//'_spec.dat')
!	do i=1,mt
!           freal = df*(i-1)         
!           o     = pi2*freal - aj*oi
!	   write(21,'(8e20.10)')real(o), imag(o), real(rrr(i)), imag(rrr(i)), 
!     &				real(fff(i)), imag(fff(i)),  
!     &				real(zzz(i)), imag(zzz(i))
!	end do
!	close(21)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ifft
	call fft(rrr,m,-1)
	call fft(fff,m,-1)
	call fft(zzz,m,-1)

!        do i=1,10
!            call fft(integt(:,i),m,-1)
!        end do
!
!        rrr=integt(:,3)*fpsv2+integt(:,4)*fpsv1+integt(:,5)*fpsv01
!     &      +integt(:,6)*fpsv02
!        fff=integt(:,1)*fsh1+integt(:,2)*fsh2
!        zzz=integt(:,7)*fpsv2+integt(:,8)*fpsv1+integt(:,9)*fpsv01
!     &      +integt(:,10)*fpsv02
!c      post-FFT process:
	print *, trim(Outname)
        open(11,file=trim(Outname)//'_1.dat')
        open(12,file=trim(Outname)//'_2.dat')
        open(13,file=trim(Outname)//'_3.dat')
        open(14,file='g_'//trim(Outname)//'.txt')
        do i=1,mt
           t=(i-1)*dt
		urr(i) = real(rrr(i))*exp(oi*t)
                uff(i) = real(fff(i))*exp(oi*t)
                uzz(i) = real(zzz(i))*exp(oi*t)
           write (11, '(2e20.10)') t, real(rrr(i))*exp(oi*t)
           write (12, '(2e20.10)') t, real(fff(i))*exp(oi*t)
           write (13, '(2e20.10)') t, real(zzz(i))*exp(oi*t)
           write (14, '(11e20.10)') t,
     &                        (real(integt(i,j))*exp(oi*t),j=1,10)
        end do
     
        close (11)
        close (12)
        close (13)
        close (14)

	call directarrival(nly,z(0:nly),vp0,vs0,r0,z0,zs,tp,ts)

        fname = trim(Outname)//'.R.sac'
        call write2sac(fname,dt,mt,urr,tp,ts)
        fname = trim(Outname)//'.T.sac'
        call write2sac(fname,dt,mt,uff,tp,ts)
        fname = trim(Outname)//'.Z.sac'
        call write2sac(fname,dt,mt,uzz,tp,ts)
	
	do j=1,10
            	call fft(integt(:,j),m,-1)
                do i=1,mt
                        t=(i-1)*dt
                        tmp(i)=real(integt(i,j))*exp(oi*t)
                end do
                write(fname,'(a,a,i2.2,a)')'g_',trim(Outname),j,'.sac'
!               print *,fname
                call write2sac(fname,dt,mt,tmp,tp,ts)
        end do
        
	end   


