module raytracing

    type Velocity
        integer :: nlayer
        real*8, allocatable :: h(:),vp(:),vs(:),z(:)
    end type Velocity

    real*8, parameter :: err = 1e-5
    
    contains
    
    subroutine ray_out(vmodel,r0,sdep,wtype,fp,ftime)
    implicit none
    type(Velocity),intent(in) :: vmodel 
    real*8,intent(in) :: r0,sdep
    character,intent(in) :: wtype
    real*8,intent(out) :: fp,ftime
    
    real*8 :: p,atime,pn,antime,pr,artime
    integer :: inflag,irflag
    
    ! for direct wave
    call direct(vmodel,r0,sdep,wtype,p,atime)
    ! for refract wave
    call refract(vmodel,r0,sdep,wtype,pn,antime,inflag)
    ! for reflect wave
    !call reflect(vmodel,r0,sdep,wtype,pr,artime,irflag)
    
    !print *, 'Direct wave :' , atime,'Reflected wave :',artime,'Refracted wave :',antime
    
    if (inflag == 0) then
        fp = p
        ftime = atime
    else
        if (antime < atime) then
            fp = pn
            ftime = antime
        else
            fp = p
            ftime = atime
        end if
    end if
    
    end subroutine ray_out
    
    subroutine model_info(modeldata,vmodel)
    implicit none
    
    character(len=*),intent(in) :: modeldata
    type(Velocity),intent(out) :: vmodel
    
    character*80 :: buffer
    integer :: nlayer,i,status
    real*8 :: rho
    
    open(10,file=modeldata,status='old',iostat=status)
    if (status /= 0) then
        print *, 'Unable to open file :', modeldata
        close(10)
        return
    end if
    
    nlayer = 0
    do while(.true.)
        read(10,fmt='(a80)',iostat=status) buffer
        if(status /= 0) exit
        nlayer = nlayer+1
    end do
    close(10)
    
    vmodel%nlayer = nlayer
    allocate(vmodel%z(nlayer))
    allocate(vmodel%h(nlayer))
    allocate(vmodel%vp(nlayer))
    allocate(vmodel%vs(nlayer))
    
    open(10,file=modeldata)
    do i=1,nlayer
        read(10,*) vmodel%z(i),vmodel%vs(i),vmodel%vp(i),rho
    end do
    close(10)
    
    vmodel%h = 0.0
    do i=1,nlayer-1
        vmodel%h(i) = vmodel%z(i+1)-vmodel%z(i)
    end do
    
    return
    end subroutine model_info
    
    subroutine direct(vmodel,theta,zs,phase,p,atime)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    real*8,intent(in) :: theta,zs
    character,intent(in) :: phase
    real*8,intent(out) :: p,atime
    
    real*8,allocatable :: z(:),h(:),v(:),the(:,:),ttime(:)
    integer :: s,i
    
    s = sourcelayer(vmodel,zs)
    allocate(z(vmodel%nlayer))
    allocate(h(vmodel%nlayer))
    allocate(v(vmodel%nlayer))
    z = vmodel%z
    h = vmodel%h
    if (phase == 'P') then
        v = vmodel%vp
    else
        v = vmodel%vs
    end if
    
    ! for source on the ground
    if ((s == 0) .and. (zs == 0.0)) then
        p = 1/0.0
        atime = theta/v(1)
        return
    end if
    
    call layered_2p_d(vmodel,theta,zs,phase,p)
    
    allocate(the(s+1,2))
    allocate(ttime(s+1))
    the = 0.0
    ttime = 0.0
    
    if (s == 1) then
        the(1,1) = theta
        the(1,2) = theta-p*v(1)*(zs-z(s))/sqrt(1-(p*v(1))**2)
        ttime(1) = sqrt((the(1,2)-the(1,1))**2+(zs-z(s))**2)/v(1)
    else
        the(1,1) = theta
        the(1,2) = theta-p*v(1)*h(1)/sqrt(1-(p*v(1))**2)
        ttime(1) = sqrt((the(1,2)-the(1,1))**2+h(1)*h(1))/v(1)
        do i=2,s-1
            the(i,1) = the(i-1,2)
            the(i,2) = the(i,1)-p*h(i)*v(i)/sqrt(1-(p*v(i))**2)
            ttime(i) = sqrt((the(i,2)-the(i,1))**2+h(i)*h(i))/v(i)
        end do
        the(s,1) = the(s-1,2)
        the(s,2) = the(s,1)-p*(zs-z(s))*v(s)/sqrt(1-(p*v(s))**2)
        ttime(s) = sqrt((the(s,2)-the(s,1))**2+(zs-z(s))**2)/v(s)
    end if
    
    atime = sum(ttime)
    
    deallocate(z)
    deallocate(h)
    deallocate(v)
    deallocate(the)
    deallocate(ttime)
    
    return    
    end subroutine direct
    
    subroutine reflect(vmodel,theta,zs,phase,pr,artime,iflag)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    real*8,intent(in) :: theta,zs
    character,intent(in) :: phase
    real*8,intent(out) :: pr,artime
    integer,intent(out) :: iflag
    
    real*8,allocatable :: z(:),h(:),v(:),the(:,:),ttime(:)
    real*8,allocatable :: p(:),atime(:)
    integer :: i,j,k,m,s,n,ns,counter,np
    
    iflag = 0
    s = sourcelayer(vmodel,zs)
    allocate(z(vmodel%nlayer))
    allocate(h(vmodel%nlayer))
    allocate(v(vmodel%nlayer))
    z = vmodel%z
    h = vmodel%h
    
    n = vmodel%nlayer
    
    if ((n == 1) .or. (s == n)) then
        print *, 'There is no reflected ',phase,' wave!'
        return
    end if
    
    if (phase == 'P') then
        v = vmodel%vp
    else
        v = vmodel%vs
    end if
    
    allocate(p(n-s))
    np = n-s
    p = 0.0
    do i=s,n-1
        call layered_2p_r(vmodel,i,theta,zs,phase,p(i-s+1))
    end do
    
    if ((n-s == 1) .and. p(1) == 0.0) then
        print *, 'There are no reflected ',phase,' waves for this source location'
        return
    end if
    
    allocate(atime(np))
    atime = 0.0
    do i=1,np
        m = s+i-1
        allocate(the(2*m-s+1,2))
        allocate(ttime(2*m-s+1))
        the = 0.0
        ttime = 0.0
        the(1,1) = theta
        the(1,2) = theta-p(i)*v(1)*h(1)/sqrt(1-(p(i)*v(1))**2)
        ttime(1) = sqrt((the(1,2)-the(1,1))**2+h(1)*h(1))/v(1)
        do k=2,m
            the(k,1) = the(k-1,2)
            the(k,2) = the(k,1)-p(i)*h(k)*v(k)/sqrt(1-(p(i)*v(k))**2)
            ttime(k) = sqrt((the(k,2)-the(k,1))**2+h(k)*h(k))/v(k)
        end do
        j = m+1
        do k=m,s+1,-1
            the(j,1) = the(j-1,2)
            the(j,2) = the(j,1)-p(i)*h(k)*v(k)/sqrt(1-(p(i)*v(k))**2)
            ttime(j) = sqrt((the(j,2)-the(j,1))**2+h(k)*h(k))/v(k)
            j = j+1
        end do
        the(2*m-s+1,1) = the(2*m-s,2)
        the(2*m-s+1,2) = the(2*m-s,2)-p(i)*(z(s+1)-zs)*v(s)/sqrt(1-(p(i)*v(s))**2)
        ttime(2*m-s+1) = sqrt((the(2*m-s+1,2)-the(2*m-s+1,1))**2+(z(s+1)-zs)**2)/v(s)
        atime(i) = sum(ttime)
        
        deallocate(the)
        deallocate(ttime)
    end do
    
    pr = 0.0
    artime = maxval(atime)+1.0
    do i=1,np
        if (artime > atime(i)) then
            pr = p(i)
            artime = atime(i)
            iflag = 1
        end if
    end do
    deallocate(z)
    deallocate(h)
    deallocate(v)
    deallocate(p)
    deallocate(atime)
    
    return
    end subroutine reflect
    
    subroutine refract(vmodel,theta,zs,phase,pn,antime,iflag)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    real*8,intent(in) :: theta,zs
    character,intent(in) :: phase
    real*8,intent(out) :: pn,antime
    integer,intent(out) :: iflag
    
    real*8,allocatable :: z(:),h(:),v(:),the(:,:),ttime(:)
    real*8,allocatable :: p(:),atime(:)
    integer,allocatable :: lidx(:),flag(:)
    integer :: i,j,k,m,s,n,ns,counter,np
    real*8 :: l_refract

    iflag = 0
    s = sourcelayer(vmodel,zs)
    allocate(z(vmodel%nlayer))
    allocate(h(vmodel%nlayer))
    allocate(v(vmodel%nlayer))
    z = vmodel%z
    h = vmodel%h
    
    n = vmodel%nlayer
    ns = s
    
    if (s == n) then
        print *, 'There is no refracted ',phase,' wave!'
        return
    end if
    
    if (phase == 'P') then
        v = vmodel%vp
    else
        v = vmodel%vs
    end if
    
    counter = 0
    do k=s,n-1
        if (v(k) < v(k+1)) counter = counter+1
    end do
    if (counter == 0) then
        print *, 'There are no refracted',phase,'waves for this source location'
        return
    else
            
        allocate(p(counter))
        allocate(lidx(counter))
        counter = 0
        do k=s,n-1
            if (v(k) < v(k+1)) then
                counter = counter+1
                p(counter) = 1/v(k+1)
                lidx(counter) = k
            end if
        end do
        np = counter
    end if
     
    allocate(atime(np))
    allocate(flag(np))
    atime = 0.0
    flag = 0
    do i=1,np
        m = lidx(i)
        allocate(the(2*m-ns+2,2))
        allocate(ttime(2*m-ns+2))
        the(1,1) = theta
        the(1,2) = theta-p(i)*v(1)*h(1)/sqrt(1-(p(i)*v(1))**2)
        ttime(1) = sqrt((the(1,2)-the(1,1))**2+h(1)*h(1))/v(1)
        do k=2,m
            the(k,1) = the(k-1,2)
            the(k,2) = the(k,1)-p(i)*h(k)*v(k)/sqrt(1-(p(i)*v(k))**2)
            ttime(k) = sqrt((the(k,2)-the(k,1))**2+h(k)*h(k))/v(k)
            ! for given p, refract wave may cannot come back to the source
            ! point. if the(k,2) is less than 0, it indicates that the
            ! raytracing fails.
            if (the(k,2) <= 0) then
                flag(i) = 1
                exit
            end if
        end do
        if (flag(i) == 1) then
            deallocate(the)
            deallocate(ttime)
            cycle
        end if
        j = m+1
        do k=m,s+1,-1
            the(j,1) = the(j-1,2)
            the(j,2) = the(j,1)-p(i)*h(k)*v(k)/sqrt(1-(p(i)*v(k))**2)
            ttime(j) = sqrt((the(j,2)-the(j,1))**2+h(k)*h(k))/v(k)
            
            if (the(j,2) <= 0) then
                flag(i) = 1
                exit
            end if
            
            j = j+1
        end do
        if (flag(i) == 1) then
            deallocate(the)
            deallocate(ttime)
            cycle
        end if
        the(2*m-s+1,1) = the(2*m-s,2)
        the(2*m-s+1,2) = the(2*m-s+1,1)-p(i)*(z(s+1)-zs)*v(s)/sqrt(1-(p(i)*v(s))**2) 
        ttime(2*m-s+1) = sqrt((the(2*m-s+1,2)-the(2*m-s+1,1))**2+(z(s+1)-zs)**2)/v(s)
        if (the(2*m-s+1,2) <= 0) then
            flag(i) = 1
            deallocate(the)
            deallocate(ttime)
            cycle
        end if
        l_refract = the(2*m-s+1,2)
        do k=2*m-s+2,m+1,-1
            the(k,:) = the(k-1,:)-l_refract
            ttime(k) = ttime(k-1)
        end do
        the(m+1,1) = the(m,2)
        the(m+1,2) = the(m+2,1)
        ttime(m+1) = l_refract*p(i)
        atime(i) = sum(ttime)
        
        deallocate(the)
        deallocate(ttime)
    end do
    
    pn = 0.0
    antime = maxval(atime)+1.0
    do i=1,np
        if ((flag(i) == 0) .and. (antime >= atime(i))) then
            pn = p(i)
            antime = atime(i)
            iflag = 1
        end if
    end do
    deallocate(z)
    deallocate(h)
    deallocate(v)
    deallocate(p)
    deallocate(lidx)
    deallocate(atime)
    deallocate(flag)
    
    return
    end subroutine refract
    
    function sourcelayer(vmodel,zs)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    real*8,intent(in) :: zs
    
    integer :: sourcelayer
    
    integer :: i,idx
    real*8,allocatable :: z(:)
    
    allocate(z(vmodel%nlayer))
    z = vmodel.z
    
    if (zs == 0) then
        idx = 1
    end if
    
    do i=2,vmodel%nlayer
        if ((zs > z(i-1)) .and. (zs <= z(i))) then
            idx = i-1
            exit
        end if
    end do
    
    if (zs > z(vmodel%nlayer)) idx = vmodel%nlayer
    
    sourcelayer = idx
    
    deallocate(z)
    return
    end function sourcelayer

    subroutine reflect_layered(vmodel,theta,zs,phase,p)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    real*8,intent(in) :: theta,zs
    character,intent(in) :: phase
    
    real*8,allocatable :: p(:)
    
    integer :: n,s,i
    
    s = sourcelayer(vmodel,zs)
    n = vmodel%nlayer
    
    if (phase == 'P') then
        ! for P wave
        if ((n == 1) .or. (s == n)) then
            print *, 'There are no reflected P waves for this source location'
            p = 0.0
            return
        end if
            
        allocate(p(n-s))
        p = 0.0
        do i=s,n-1
            call layered_2p_r(vmodel,i,theta,zs,phase,p(i))
        end do
        
    else
        ! for S wave
        if ((n == 1) .or. (s == n-1)) then
            print *, 'There are no reflected S waves for this source location'
            p = 0.0
            return
        end if
            
        allocate(p(n-s))
        p = 0.0
        do i=s,n-1
            call layered_2p_r(vmodel,i,theta,zs,phase,p(i))
        end do
                
    end if
        
    return
    end subroutine reflect_layered
    
    ! for direct wave, n=s; for reflected wave, z(n) is the reflected plane
    ! string is either 'd' or 'r'
    ! for direct P/S wave  
    subroutine layered_2p_d(vmodel,theta,zs,phase,p)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    real*8,intent(in) :: theta,zs
    character,intent(in) :: phase
    
    real*8,intent(out) :: p
    
    integer :: k,s,n,M,counter
    real*8,allocatable :: h(:),v(:),z(:),hh(:)
    real*8 :: q,err1
    
    s = sourcelayer(vmodel,zs)
    
    allocate(z(vmodel%nlayer))
    allocate(h(vmodel%nlayer))
    allocate(v(vmodel%nlayer))
    
    if (phase == 'P') then
        v = vmodel.vp
    else
        v = vmodel.vs
    end if
    
    h = vmodel.h
    z = vmodel.z
    n = s
    
    allocate(hh(n))
    hh(1:s-1) = h(1:s-1)
    hh(s) = zs-z(s)
    
    M = 1
    do k=1,n
        if (v(k) > v(M)) M = k
    end do
    
    q = initial(n,M,v,hh,theta)
    err1 = abs(F(q,n,M,v,hh,theta))
    counter = 0
    do while (err1 > err)
        q = q-F(q,n,M,v,hh,theta)/F1(q,n,M,v,hh)
        err1 = abs(F(q,n,M,v,hh,theta))
        counter = counter+1
        if (counter > 200) print *, 'for direct wave, q is unstable'
    end do
    
    p = q/(v(M)*sqrt(1+q*q))
    
    deallocate(z)
    deallocate(h)
    deallocate(v)
    deallocate(hh)
    return
    end subroutine layered_2p_d
    
    ! for reflect P/S wave
    subroutine layered_2p_r(vmodel,n,theta,zs,phase,p)
    implicit none
    
    type(Velocity),intent(in) :: vmodel
    integer,intent(in) :: n
    real*8,intent(in) :: theta,zs
    character,intent(in) :: phase
    
    real*8,intent(out) :: p
    
    integer :: k,s,M,counter
    real*8,allocatable :: h(:),v(:),z(:),hh(:)
    real*8 :: q,err1
    
    allocate(z(vmodel%nlayer))
    allocate(h(vmodel%nlayer))
    allocate(v(vmodel%nlayer))
    
    s = sourcelayer(vmodel,zs)
    if (phase == 'P') then
        v = vmodel.vp
    else
        v = vmodel.vs
    end if
    
    h = vmodel.h
    z = vmodel.z
    
    allocate(hh(n))
    hh(1:s-1) = h(1:s-1)
    hh(s) = z(s+1)-zs+h(s)
    hh(s+1:n) = 2*h(s+1:n)
    
    M = 1
    do k=1,n
        if(v(k) > v(M)) M = k
    end do
    
    q = initial(n,M,v,hh,theta)
    err1 = abs(F(q,n,M,v,hh,theta))
    counter = 0
    do while (err1 > err)
        q = q-F(q,n,M,v,hh,theta)/F1(q,n,M,v,hh)
        err1 = abs(F(q,n,M,v,hh,theta))
        counter = counter+1
        if (counter > 200) print *, 'for reflect wave, q is unstable'
    end do
    
    p = q/(v(M)*sqrt(1.0+q*q))
    
    deallocate(z)
    deallocate(h)
    deallocate(v)
    deallocate(hh)
    return
    end subroutine layered_2p_r
    
    function F(q,n,M,v,hh,theta) 
    implicit none
    
    integer,intent(in) :: n,M
    real*8,intent(in) :: q,v(n),hh(n),theta
    
    real*8 :: F
    integer :: k
    
    F = 0.0
    do k=1,n
        F = F+hh(k)*v(k)/sqrt(v(M)*v(M)+q*q*(v(M)*v(M)-v(k)*v(k)))
    end do
    F = q*F-theta
    
    return
    end function F
    
    function F1(q,n,M,v,hh) 
    implicit none
    
    integer,intent(in) :: n,M
    real*8,intent(in) :: q,v(n),hh(n)
    
    real*8 :: F1
    integer :: k
    
    F1 = 0.0
    do k=1,n
        F1 = F1+hh(k)*v(k)/(v(M)*v(M)+q*q*(v(M)*v(M)-v(k)*v(k)))**1.5
    end do
    F1 = F1*v(M)*v(M)
    
    return
    end function F1
    
    function initial(n,M,v,hh,theta)
    implicit none
    
    integer,intent(in) :: n,M
    real*8,intent(in) :: v(n),hh(n),theta
    
    real*8 :: initial
    
    integer :: k
    real*8 :: a,b,etak,thetac,q0
    
    a = 0.0
    b = 0.0
    do k=1,n
        etak = v(k)/v(M)
        a = a+etak*hh(k)/hh(M)
        
        if (k /= M) b = b+etak*hh(k)/(1-etak*etak)
    end do
    
    thetac = a*b/(a-1.0)
    if (theta < thetac) then
        q0 = theta/a
    else
        q0 = theta-b
    end if
    
    initial = q0
    
    return
    end function initial
    
    
end module raytracing
    
