subroutine reftrace(alfa,zlow,zup,ips, time,dist,hmax)


! zstat : starting level
! alfa : dipping angle at the starting level
! zzt : finishing level. It is always above that zstart
! izgib : marker. When 0, the ray goes

real zzz(300),rrr(300),vvv(300)

common/refmod/nrefmod,zref(600),vref(600,2)
common/pi/pi,per
rz=6378.

!write(*,*)alfa,zstat,zzt,izgib

!zup=-5
!zlow=-1
!ips=1
!alfa=60

dist=-999.
time=-999.

rup=rz-zup
rlow=rz-zlow

vup=vrefmod(zup,ips)
vlow=vrefmod(zlow,ips)
!write(*,*)' zzt=',zzt,' vzt=',vzt,' ips=',ips


nref=nrefmod+1
do i=1,nref-1
	rrr(i)=rz-zref(i)
	vvv(i)=vref(i,ips)
end do
rrr(nref)=0
vvv(nref)=15.

!do i=1,nref
!    write(*,*)i,rz-rrr(i),vvv(i)
!end do
!stop

if(abs(alfa-180).lt.1.e-9) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Tracing vertical ray from rlow to rup
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    time=0
    hmax=zlow
    dist=0
    do il=1,nref-1
        v1=vvv(il)
        v2=vvv(il+1)
        r1=rrr(il)
        r2=rrr(il+1)
        !write(*,*)' ini: r1=',rz-r1,' r2=',rz-r2
        if(r2.ge.rup-1.e-9) cycle
        if(r1.le.rlow) exit
        if((r1-rup)*(r2-rup).lt.0.)then
            v1=vup
            r1=rup
        end if
        if((r1-rlow)*(r2-rlow).le.0.)then
            v2=vlow
            r2=rlow
        end if
        !write(*,*)' aft: r1=',rz-r1,' r2=',rz-r2
        !write(*,*)' aft: v1=',v1,' v2=',v2
        p2=r2/v2
        if(r2.eq.r1.and.p2.gt.p) cycle
        if(r2.eq.r1.and.p2.le.p) then
            hmax=rz-r2
            exit
        end if

        aaa = (v2-v1) / (r2-r1)
        bbb = v1 - aaa*r1
        dt = (alog(aaa*r1+bbb) - alog(aaa*r2+bbb)) / aaa
        time=time+dt

    end do
    return
end if

sina=sin(alfa*per)
p=rlow*sina/vlow

!write(*,*)' p=',p,' vzt=',vzt


!write(*,*)' rup=',rup,' rlow=',rlow

ep1=0.
time1=0.

k_1_2 = 2

if(abs(alfa-90).lt.1.e-9) then ! The case when ray starts horizontally
    rmax=rlow
    hmax=rz-rmax
    p=rlow/vlow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute the maximal depth of the ray
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do ilay=1,nref-1
!	v1=vvv(ilay)
!	v2=vvv(ilay+1)
!	r1=rrr(ilay)
!	r2=rrr(ilay+1)
!	p2=r2/v2
!	if(p2.gt.p) cycle
!	if(r2.eq.r1.and.p2.lt.p) exit
!	aaa=(v1-v2)/(r1-r2)
!	bbb=(v2*r1-v1*r2)/(r1-r2)
!	rmax=p*bbb/(1.-p*aaa)
!	if((r1-rmax)*(r2-rmax).le.0.)goto 53
!end do
!ilay=ilay-1
!53 continue
!hmax=rz-rmax
!vmax=vrefmod(hmax,ips)
!
!write(*,*)' hmax=',hmax,' vmax=',vmax





    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Tracing from rlow to the turning point
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    time=0
    ep=0
    do il=1,nref-1
	    v1=vvv(il)
	    v2=vvv(il+1)
	    r1=rrr(il)
	    r2=rrr(il+1)
    !write(*,*)' ini: r1=',rz-r1,' r2=',rz-r2
	    if(r2.ge.rup-1.e-9) cycle
	    if(r1.le.rmax) exit
	    if((r1-rup)*(r2-rup).lt.0.)then
		    v1=vup
		    r1=rup
	    end if
	    if(r2.eq.r1) cycle
            !write(*,*)' aft: r1=',rz-r1,' r2=',rz-r2
	    call ray_lay(r1,r2,v1,v2,p,dt,dep)
	    !write(*,*)' v1=',v1,' v2=',v2,' p=',p
	    !write(*,*)' dt=',dt,' dep=',dep*rz
	    time=time+dt
	    ep=ep+dep
    !write(*,*)' dep=',dep,' dt=',dt
    !write(*,*)
    end do
    dist=ep*rz
!write(*,*)' ep2=',ep*rz,' time2=',time

    return
end if


if(abs(rup-rlow).lt.1.e-9) then ! The case when source and receiver are at the same level
    k_1_2 = 2
    goto 883
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tracing from rlow to rup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do il=1,nref-1
	v1=vvv(il)
	v2=vvv(il+1)
	r1=rrr(il)
	r2=rrr(il+1)
!write(*,*)' ini: r1=',rz-r1,' r2=',rz-r2
	if(r2.ge.rup-1.e-9) cycle
	if(r1.le.rlow) exit
	if((r1-rup)*(r2-rup).lt.0.)then
		v1=vup
		r1=rup
	end if
	if((r1-rlow)*(r2-rlow).le.0.)then
		v2=vlow
		r2=rlow
	end if
!write(*,*)' aft: r1=',rz-r1,' r2=',rz-r2
!write(*,*)' aft: v1=',v1,' v2=',v2
	p2=r2/v2
	if(r2.eq.r1.and.p2.gt.p) cycle
	if(r2.eq.r1.and.p2.le.p) then
		hmax=rz-r2
		exit
	end if
	call ray_lay(r1,r2,v1,v2,p,dt,dep)
!write(*,*)' z1=',rz-r1,' z2=',rz-r2,dt,dep*rz
	!if(il.eq.16)write(*,*)' v1=',v1,' v2=',v2,' p=',p
	!write(*,*)il,' dt=',dt,' dep=',dep,' dz=',r2-r1
	!if(il.eq.32) pause
	time1=time1+dt
	ep1=ep1+dep
	if(p2.le.p) then
		aaa=(v1-v2)/(r1-r2)
		bbb=(v2*r1-v1*r2)/(r1-r2)
		hmax=rz-p*bbb/(1.-p*aaa)
		!write(*,*)' h1=',rz-r1,' h2=',rz-r2
		!write(*,*)' v1=',v1,' v2=',v2
		!write(*,*)' hmax=',hmax
		return
	end if
	!write(*,*)' il=',il,' h1=',rz-r1,' h2=',rz-r2
	!if(il.eq.16)write(*,*)il,v1,v2,ep1*rz,dep*rz
end do

!write(*,*)' ep1=',ep1*rz,' time1=',time1

883 continue

ep2=0.
time2=0.
hmax=zlow

if(alfa.gt.90+1.e-9)goto 334

!write(*,*)' rup=',rup,' rlw=',rlow
!write(*,*)' vup=',vup,' vlw=',vlow


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the maximal depth of the ray
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ilay=1,nref-1
	v1=vvv(ilay)
	v2=vvv(ilay+1)
	r1=rrr(ilay)
	r2=rrr(ilay+1)
	p2=r2/v2
!write(*,*)' z1=',rz-r1,' z2=',rz-r2
!write(*,*)' v1=',v1,' v2=',v2
!write(*,*)' p=',p,' p2=',p2
	if(p2.gt.p) cycle
	if(r2.eq.r1.and.p2.lt.p) exit
	aaa=(v1-v2)/(r1-r2)
	bbb=(v2*r1-v1*r2)/(r1-r2)
	rmax=p*bbb/(1.-p*aaa)
        if(rmax.lt.r2) rmax=r2
    !write(*,*)' zmax=',rz-rmax
	!if((r1-rmax)*(r2-rmax).le.0.)goto 3
	goto 3
end do
ilay=ilay-1
3 continue
hmax=rz-rmax
vmax=vrefmod(hmax,ips)

!write(*,*)' hmax=',hmax,' vmax=',vmax



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tracing from rlow to the turning point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
884	continue

do il=1,nref-1
	v1=vvv(il)
	v2=vvv(il+1)
	r1=rrr(il)
	r2=rrr(il+1)
!write(*,*)' ini: r1=',rz-r1,' r2=',rz-r2
	if(r2.ge.rlow-1.e-9) cycle
	if(r1.le.rmax) exit
	if((r1-rlow)*(r2-rlow).lt.0.)then
		v1=vlow
		r1=rlow
	end if
	if(r2.eq.r1) cycle
!write(*,*)' aft: r1=',rz-r1,' r2=',rz-r2
	call ray_lay(r1,r2,v1,v2,p,dt,dep)
	!write(*,*)' v1=',v1,' v2=',v2,' p=',p
	!write(*,*)' dt=',dt,' dep=',dep*rz
	time2=time2+dt
	ep2=ep2+dep
!write(*,*)' dep=',dep,' dt=',dt
!write(*,*)
end do

!write(*,*)' ep2=',ep2*rz,' time2=',time2

334 continue


time = time2 * 2 + time1
ep = ep2 * 2 + ep1
dist=ep*rz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*)' time=',time,' dist=',dist
return
end


