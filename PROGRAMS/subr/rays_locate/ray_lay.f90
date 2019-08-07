subroutine ray_lay(r2,r1,v2,v1,p,time,dist)
!     This subroutine computes transit times and delta increase
!     of a ray travelling from radius R2 down to R1 (or to the
!     turning point Rt if R1&lt;Rt&lt;R2), for linear velocity behaviour
!     from V2 (at R2) to V1. Slowness is P. Output: time T (sec)
!     and distance D (rad). Velocity in km/sec, radii in km.
!     Equations from G.Nolet, Linearized inversion of (teleseismic)
!     data. In: R.Cassinis,'The inverse problem in geophysical
!     interpretation',Plenum Press,NY,1981.
common/pi/pi,per
rz=6378

rr=(r2-r1)/r1
vv=(v2-v1)/v1
!kod_sing=0
!if(abs(vv/rr)-1..lt.0.01) then
!	kod_sing=1
!	rzz=rz+1000.
!	p = p * rzz / rz
!	r1=r1+1000.
!	r2=r2+1000.
!	write(*,*)' singularity'
!end if
!write(*,*)' kod_sing=',kod_sing


time=0
dist=0
if(abs(v1*v2).lt.0.000001) return

if(abs(p).lt.0.000001) then
	if(abs(v2-v1).gt.0.000001) then
		time=(log(v2)-log(v1))/dvdr
	else
		time=(r2-r1)/v1
	end if
	return
end if


!vr1=v1/r1
!vr2=v2/r2
!dif_pr=100*(vr2-vr1)/vr1
!write(*,*)' v1/r1=',v1/r1,' v2/r2=',v2/r2,' dif_pr=',dif_pr
!if(abs(dif_pr).lt.0.0001) then
!	sint=p*v1/r1
!	cost=sqrt(1-sint*sint)
!	ds=abs(r2-r1)/cost
!	v=(v2+v1)/2
!	r=(r2+r1)/2
!	time=ds/v
!	dist=ds*sint/r
!	return
!end if


zav=rz-(r2+r1)/2.

!if((zav-150)*(zav-155).lt.0) write(*,*)' z1=',rz-r1,' z2=',rz-r2
!if((zav-150)*(zav-155).lt.0) write(*,*)' v1=',v1,' v2=',v2
!if((zav-150)*(zav-155).lt.0) write(*,*)' v1/r1=',v1/r1,' v2/r2=',v2/r2

!write(*,*)
!write(*,*)' v1/r1=',v1/r1,' v2/r2=',v2/r2,' diff=',diff

diff=100*(v1/r1-v2/r2)/(v1/r1)
if(diff.lt.0.0001) then
	sint=p*v1/r1
	cost=sqrt(1-sint*sint)
	dz=abs(r2-r1)
	ds=dz/cost
	vav=(v1+v2)/2
	time=ds/vav
	dh=ds*sint
	dist=dh/rz
	return
end if



dvdr=(v2-v1)/(r2-r1)
ap = dvdr*p

if(r2.gt.0.0000001) then
	sint2=p*v2/r2
else
	sint2=1.
end if
if(sint2.gt.1.) return

if(r1.gt.0.0000001) then
	sint1=p*v1/r1
else
	sint1=1.
end if
if(sint1.gt.1.) sint1=1


!write(*,*)' sint1=',sint1,' sint2=',sint2

teta1=asin(sint1)
teta2=asin(sint2)

tgt1=TAN(teta1 / 2.)
tgt2=TAN(teta2 / 2.)

if(ap.eq.0.) then
	sntet=0
	c1=1. / TAN(TETA1)
	c2=1. / TAN(TETA2)
	time=p*(c2-c1)
else
	ap2=ap*ap
	sqap=SQRT(abs(1.-ap2))
	if(ap2.eq.1.) then
		sntet1=TAN(.785398132 + ap * .5 * teta1)
		sntet2=TAN(.785398132 + ap * .5 * teta2)
		!write(*,*)' metka2'
	else if(ap2.lt.1.) then
		alfa=acos(-ap)

		sntet1=-.5*log((1.+SIN(alfa+teta1))/(1.-SIN(alfa-teta1)))/sqap
		sntet2=-.5*log((1.+SIN(alfa+teta2))/(1.-SIN(alfa-teta2)))/sqap
		!write(*,*)' teta1=',teta1,' teta2=',teta2,' alfa=',alfa
	else if(ap2.gt.1) then
		aptg1 = -ap * tgt1 + 1.
		aptg2 = -ap * tgt2 + 1.
		sntet1=2.*ATAN(aptg1/sqap)/sqap
		sntet2=2.*ATAN(aptg2/sqap)/sqap
		!write(*,*)' metka222'
	end if
	time = (log(tgt2) - log(tgt1) - sntet2 + sntet1) / dvdr
	dist = teta1-teta2-ap*(sntet2-sntet1)
	!write(*,*)' metka11'
end if

!if(kod_sing.eq.1) then
!	dist=dist*rzz/rz
!	p = p * rz / rzz
!	r1=r1-1000.
!	r2=r2-1000.
!end if


return
end