!ccccccccccccccc Discrete Wavenumber Integration Method cccccccccccccccccccc

      subroutine DWIM(o, nf, kmin, kmax, sums)
      
      use comvar
      include     "green_com.inc"
      
      integer    :: nf, nk
      real*8     :: kmin, kmax, restdk
      complex*16 :: integ(10), sums(10), integ1(10), integ2(10)
      
      integer :: counter
      integer :: flagsum
      
! use trapezoidal integration and assume that the first step is strat from k=0
      kn = kmin
      call integrand_calc(kn, o, integ1)     
      if (kmin.eq.0.0) then
            sums = integ1*dk/2d0
      else
            sums = integ1*dk
      end if          
      
      nk = floor(kmax/dk)+1
      do i=1,nk
            kn = kmin+i*dk
            
            call integrand_calc(kn, o, integ2)
            
            sums = sums + integ2*dk
!            integ1 = integ2

      end do      

!      restdk = dk-(kn-kmax)
!      kn = kmax
!      call integrand_calc(kn, o, integ2)
!      sums = sums + (integ1/2d0+integ2)*restdk
                
! because integrand oscillates and slowly attenuate, Peak Trough Averaging Method 
! has been used to accelerate the attenuation and improve precision of integration

      counter = 0
      flagrp = 1        !global array
      flagip = 1        !global array
      num = 0           !global array
      flag = 1          !global array
      pre = 0d0         !global array
      tre = 0d0         !global array
      pre = 0d0         !global array
      tim = 0d0         !global array
      rp = 0d0          !global array
      ip = 0d0          !global array
      
!           write (*, '(a2,10e20.10)')'in',(real(sums(j)),j=1,10)
      do 
            counter = counter+1
            kn = kn+dk
! compute integrand            
            call integrand_calc(kn, o, integ)
! go one step ahead for integration            
            sums = sums + integ*dk

! start Peak Trough Averaging Method 
            flagsum = 0             !global array
            do i=1,10
                  call ptam(i, nf, kn, sums(i), counter, o)
                  flagsum = flagsum + flagrp(i) + flagip(i)
            end do
           
            if(flagsum == 0) exit
!            kn = kn+dk
!            write(*,*) counter,nf,kmax,kn
      end do
      
      sums = sumreal + aj*sumimag
!           write (*, '(10e20.10)')((sumreal(j)),j=1,10)
      end subroutine DWIM
