subroutine halftrace(p,zstart,zbot,ips,dzstep,zmaxst, time,dist,hmax)


real zzz(8000),rrr(8000),vvv(8000)
real htmp(8000),vtmp(8000)				! ref. models in different zones
real hmod(8000),vmod(8000)				! ref. models in different zones
real hmod0(600),vmodp0(600),vmods0(600) ! ref. models in different zones
real dray(5000),hray(5000),tray(5000)


common/refmod/nrefmod0,hmod0,vmodp0,vmods0
common/ray/npath,dray,hray,tray


rz=6371.
one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

rbot=rz-zbot
hmax=zbot

nrefmod=nrefmod0

do i=1,nrefmod
	hmod(i)=hmod0(i)
	vmod(i)=vmodp0(i)
	if(ips.eq.2)vmod(i)=vmods0(i)
end do

il=1
htmp(il)=hmod(1)
vtmp(il)=vmod(1)
!write(*,*)' dzstep=',dzstep,' zmaxst=',zmaxst
!pause

do i=2,nrefmod
	h1=hmod(i-1)
	h2=hmod(i)
	v1=vmod(i-1)
	v2=vmod(i)
!	write(*,*)' h1=',h1,' h2=',h2
	if(h2-h1.gt.dzstep) then
		nzst=(h2-h1)/dzstep+1
		dzst=(h2-h1)/nzst
		dvst=(v2-v1)/nzst
		do ii=1,nzst
			il=il+1
			htmp(il)=h1+ii*dzst
			vtmp(il)=v1+ii*dvst
			if(h1+ii*dzst.gt.zmaxst) then
				il=il+1
				htmp(il)=h2
				vtmp(il)=v2
				goto 739
			end if
		end do
	else
		il=il+1
		htmp(il)=h2
		vtmp(il)=v2
!		write(*,*)il,h2,v2
		if(h2.gt.zmaxst) goto 739
	end if
end do

739 continue
hmod=htmp
vmod=vtmp
nrefmod=il


!write(*,*)' nrefmod=',nrefmod
!do i=1,nrefmod
!	write(*,*)hmod(i),vmod(i)
!end do
zend=hmod(nrefmod)

!write(*,*)' zstart=',zstart
vstart=vrefmod(zstart,ips)
!write(*,*)' vstart=',vstart

rstart=rz-zstart

nref=nrefmod
do i=1,nref
	rrr(i)=rz-hmod(i)
	vvv(i)=vmod(i)
end do

!do i=1,nref
	!write(*,*)6371-rrr(i),vvv(i)
!end do
!pause
!write(*,*)' ep1=',ep1,' time1=',time1

883 continue

ep2=0.
time2=0.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the maximal depth of the ray
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ilay=1,nref-1
	v1=vvv(ilay)
	v2=vvv(ilay+1)
	r1=rrr(ilay)
	r2=rrr(ilay+1)

	p2=r2/v2
	if(p2.gt.p) cycle
	if(r2.eq.r1.and.p2.le.p) then
		rmax=r2
		goto 3
	end if
	!write(*,*)6371-r1,6371-r2,v1,v2
	!write(*,*)' p=',p,' p2=',p2
	aaa=(v1-v2)/(r1-r2)
	bbb=(v2*r1-v1*r2)/(r1-r2)
	rmax=p*bbb/(1.-p*aaa)
	!write(*,*)' r1=',r1,' r2=',r2,' rmax=',rmax
	!pause
	goto 3
end do
!ilay=ilay-1
rmax=rz-zbot
3 continue
if(rmax.lt.rz-zmaxst)rmax=rz-zmaxst
hmax=6371.-rmax
!write(*,*)' hmax=',hmax,' zstart=',zstart

if(hmax.lt.dzstep) then
	npath=0
	return
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tracing from rlow to rmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
884	continue

npath=0
do il=1,nref-1
	v1=vvv(il)
	v2=vvv(il+1)
	r1=rrr(il)
	r2=rrr(il+1)
!write(*,*)' ini: r1=',r1,' r2=',r2
	if(r2.ge.rstart) cycle
	if(r1.le.rmax) exit
	if(r1.le.rbot) exit
	if((r1-rstart)*(r2-rstart).lt.0.)then
		v1=vstart
		r1=rstart
	end if
	if((r1-rbot)*(r2-rbot).lt.0.)then
		r2=rbot
		v2=vrefmod(zbot,ips)
	end if
	if(abs(r2-r1).lt.1.e-8) cycle
!write(*,*)' aft: r1=',r1,' r2=',r2
	call ray_lay(r1,r2,v1,v2,p,dt,dep)
!	write(*,*)
!	write(*,*)' v1=',v1,' v2=',v2,' p=',p
!	write(*,*)' r1=',r1,' r2=',r2,' dt=',dt,' dep=',dep
	time2=time2+dt
	ep2=ep2+dep
	if((r1-rmax)*(r2-rmax).lt.0.)r2=rmax

	zcur=rz-r2
	npath=npath+1
	dray(npath)=ep2*rz
	hray(npath)=zcur
	tray(npath)=time2
!	write(*,*)zcur,dep,ep2*rz,time2
	if(zcur.gt.zmaxst) goto 334

!write(*,*)' dep=',dep,' dt=',dt
!write(*,*)
end do


!write(*,*)' ep2=',ep2,' time2=',time2

334 continue

time=time2
ep=ep2
dist=ep*rz

!npath=npath+1
!dray(npath)=dist
!hray(npath)=hmax
!tray(npath)=time2

!write(*,*)' time1=',time1,' ep1=',ep1
!write(*,*)' time2=',time2,' ep2=',ep2
!write(*,*)' alfa=',alfa,' dist=',dist,' time=',time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,*)' time=',time,' dist=',dist
return
end
