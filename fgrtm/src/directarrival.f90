	subroutine directarrival(nlayer,dep,vp,vs,r0,rdep,sdep,tp,ts)
	use raytracing
	implicit none

	type(Velocity) :: vel
	integer :: nlayer
	real*8 :: dep(nlayer),vp(nlayer),vs(nlayer)
	real*8 :: rdep, sdep, r0, tp, ts
	integer :: i, iloc 
    real*8 :: fp, fs
    real*8 :: tmp

    if (rdep == sdep) then
        do i=1,nlayer
            if (sdep < dep(i)) then
                iloc = i - 1
                exit
            end if
        end do
        tp = r0/vp(iloc)
        ts = r0/vs(iloc)
    else
        if (rdep > sdep) then
            tmp = sdep
            sdep = rdep
            rdep = tmp
        end if
        
        do i=1,nlayer
            if (rdep < dep(i)) then
                iloc = i - 1
                exit
            end if
        end do
    
	    vel%nlayer = nlayer - iloc +1
	    allocate(vel%h(vel%nlayer))
        allocate(vel%z(vel%nlayer))
	    allocate(vel%vp(vel%nlayer))
	    allocate(vel%vs(vel%nlayer))

	    do i=iloc,nlayer
		    vel%vp(i-iloc+1) = vp(i)
		    vel%vs(i-iloc+1) = vs(i)
		    vel%z(i-iloc+1) = dep(i) - rdep
        end do
        vel%z(1) = 0.0

        vel%h = 0.0
        do i=1,vel%nlayer-1
            vel%h(i) = vel%z(i+1)-vel%z(i)
        end do
    !ray_out(vmodel,r0,sdep,wtype,fp,ftime)
	    call ray_out(vel,r0,sdep - rdep,'P',fp,tp)
	    call ray_out(vel,r0,sdep - rdep,'S',fs,ts)
    end if
    
	return
	end subroutine directarrival
