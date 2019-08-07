subroutine ray_xyz(xzt,yzt,zzt, xst,yst, alfa,iiips, hmax)
real dray(5000),hray(5000),tray(5000)
real dray1(5000),hray1(5000)
real dray2(5000),hray2(5000)


common /ray_path/npray,xray(5000),yray(5000),zray(5000)
common/pi/pi,per
common/ray/npath,dray,hray,tray
common/grid/zgrmax,dzlay,dsmin

dx=xzt-xst
dy=yzt-yst
dshor=sqrt(dx*dx+dy*dy)

if(dshor.lt.10.) then
	dist=sqrt(dshor*dshor + zzt*zzt)
	npp=dist/dsmin
	ds=dist/npp
	do i=1,npp+1
		xray(i)=xzt+(i-1)*ds*((xst-xzt)/dist)
		yray(i)=yzt+(i-1)*ds*((yst-yzt)/dist)
		zray(i)=zzt+(i-1)*ds*((0. -zzt)/dist)
	end do
	npray=npp+1
	return
end if


vzt=vrefmod(zzt,iiips)

sina=sin(alfa*per)
cosa=cos(alfa*per)

cosb=(dx)/dshor
sinb=(dy)/dshor


zst=0
rzt=6371.-zzt
px=rzt*sina/vzt

zbot=zgrmax
if(alfa.ge.90) zbot=zzt

!write(*,*)' alfa=',alfa
!write(*,*)' px=',px,' zst=',zst,' zbot=',zbot
!write(*,*)' dshor=',dshor

call halftrace(px,zst,zbot,iiips,dzlay,zgrmax, time,dist1,hmax)
if(hmax.gt.zbot) hmax=zbot
!do i=1,npath
	!write(*,*)' dray=',dray(i),' hray=',hray(i)
!end do

if(npath.eq.0) then
	dist=sqrt(dshor*dshor + zzt*zzt)
	npp=dist/dsmin
	ds=dist/npp
	do i=1,npp+1
		xray(i)=xzt+(i-1)*ds*((xst-xzt)/dist)
		yray(i)=yzt+(i-1)*ds*((yst-yzt)/dist)
		zray(i)=zzt+(i-1)*ds*((0. -zzt)/dist)
	end do
	npray=npp+1
	return
end if
!write(*,*)' dist=',dist1,' time=',time,' hmax=',hmax

np=1
dray1(np)=0
hray1(np)=zst

do ii=1,npath
	np=np+1
	dray1(np)=dray(ii)
	hray1(np)=hray(ii)
	!write(*,*)np,dray1(np),hray1(np)
end do
!pause
npath1=np
dist=dray(npath)


if(alfa.ge.90) then
	np=1
	dray(np)=0
	hray(np)=zzt
	do i=npath1,1,-1
		d1=dray1(i)
		h1=hray1(i)
		if(d1.ge.dshor) cycle
		if(abs(h1-zzt).lt.0.00001) cycle
		np=np+1
		dray(np)=dshor-d1
		hray(np)=h1
!				write(*,*)np,dray(np),hray(np)
	end do
	npath=np
	goto 991
end if

call halftrace(px,zzt,zbot,iiips,dzlay,zgrmax, time,dist2,hmax)

!write(*,*)' dist=',dist2,' time=',time,' hmax=',hmax

np=1
dray2(np)=0
hray2(np)=zzt

do ii=1,npath
	np=np+1
	dray2(np)=dray(ii)
	hray2(np)=hray(ii)
	!write(*,*)np,dray2(np),hray2(np)
end do

npath2=np

i2=1
h2=hray2(i2)
d2=dray2(i2)
hlast2=hray2(1)
do i1=1,npath1
	h1=hray1(i1)
	d1=dray1(i1)
	if(h1.gt.h2) then
		i2=i2+1
		h2=hray2(i2)
		d2=dray2(i2)
	end if
	if(d1+d2.gt.dshor) exit
end do
i2=i2-1
!write(*,*)' i2=',i2
np=0
do i=1,i2
	np=np+1
	dray(np)=dray2(i)
	hray(np)=hray2(i)
	!write(*,*)dray2(i),hray2(i)
end do
i1=i1-1
!write(*,*)' i1=',i1
do i=i1,1,-1
	np=np+1
	dray(np)=dshor-dray1(i)
	hray(np)=hray1(i)
	!write(*,*)dshor-dray1(i),hray1(i)
end do
!		do i=1,np
!			write(*,*)dray(i),hray(i)
!		end do
!		write(*,*)' dshor=',dshor
npath=np
!		pause

991		continue

np=1
dray2(1)=dray(1)
hray2(1)=hray(1)
do ii=2,npath
	d1=dray2(np)
	z1=hray2(np)
	d2=dray(ii)
	z2=hray(ii)
	ds=sqrt((d2-d1)*(d2-d1)+(z2-z1)*(z2-z1))	
	if(ds.lt.dsmin.and.ii.ne.npath.and.z2.lt.zgrmax) cycle
	ns=ds/dsmin	
	!write(*,*)' ds=',ds,' ns=',ns
	if(ns.eq.0)ns=1
	dz=(z2-z1)/ns	
	dd=(d2-d1)/ns
	if(z1.gt.zgrmax-10.and.z2.gt.zgrmax-10.and.ns.ne.1) then
!			write(*,*)' z1=',z1,' d1=',d1
!			write(*,*)' z2=',z2,' d2=',d2,' ns=',ns
		np=np+1
		dray2(np)=d2
		hray2(np)=z2
		cycle
	end if

	do iii=1,ns
		np=np+1
		dray2(np)=d1+iii*dd
		hray2(np)=z1+iii*dz
	end do
end do
npath2=np

!write(*,*)' npath2=',npath2
!write(*,*)' ztr:',xzt,yzt,zzt
np=0
do ii=1,npath2
	!write(*,*)ii,dray2(ii),hray2(ii)
	x=xzt-dray2(ii)*cosb
	y=yzt-dray2(ii)*sinb
	z=hray2(ii)
	np=np+1
	xray(np)=x
	yray(np)=y
	zray(np)=z
	!write(*,*)x,y,z
end do
npray=np
!write(*,*)' stat:',xst,yst,zst

return
end