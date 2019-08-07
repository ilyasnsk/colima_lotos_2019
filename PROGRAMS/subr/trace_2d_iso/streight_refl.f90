subroutine streight_refl(xzt,yzt,xst,yst,ips,	tout,xrefl,yrefl)

real xrmin(1000),yrmin(1000)
real t_segm(1000)
integer indexx(100)
real dxtmp(1000),dytmp(1000),xtmp(1000),ytmp(1000)
integer k_around(100)
real t_around(100)

common/ray_param/ds_ini,ds_segm_min,bend_min0,bend_max0

common/ray/ nodes,xray(1000),yray(1000)
common/ray_part/ npart,xpart(1000),ypart(1000),spart(1000),kod_part(1000)
common/shift/ dxpart(1000),dypart(1000)

!bend_min0=0.0001

dxmax_refl=5
dxmin_refl=0.5
nx_refl=5


!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
! Preliminary reflected ray, as two line:

if(abs(xst-xzt).lt.0.0001) then
	xrefl0=xst
else
	yzt2 = yyy_surf(xzt)*2 - yzt
	yst2 = yyy_surf(xst)*2 - yst

	!write(*,*)' yzt2=',yzt2,' yst2=',yst2

	aaa1=(yzt-yst2)/(xzt-xst)
	bbb1=(xzt*yst2-xst*yzt)/(xzt-xst)

	aaa2=(yzt2-yst)/(xzt-xst)
	bbb2=(xzt*yst-xst*yzt2)/(xzt-xst)

	xxx = (bbb2-bbb1) / (aaa1-aaa2)
	yyy1 = aaa1 * xxx + bbb1
	yyy2 = aaa2 * xxx + bbb2

	!write(*,*)' xxx=',xxx,' yyy=',yyy1,yyy2
	xrefl0=xxx
end if

yrefl0=yyy_surf(xxx)



dx_refl=dxmax_refl

k_around=0
t_around=0
nk=0
icen=0
tmin=99999

135 continue
index=0
do icur=-nx_refl,nx_refl,1

	xrefl = xrefl0 + (icur+icen) * dx_refl
	yrefl = yyy_surf(xrefl)
	if(nk.ne.0) then
		do ik=1,nk
			if(k_around(ik).eq.icur+icen) goto 134
		end do
	end if

	call streight_line(xzt,yzt,xrefl,yrefl,ips, tout1)
	call streight_line(xrefl,yrefl,xst,yst,ips, tout2)
	tout=tout1+tout2
	nk=nk+1
	k_around(nk)=icur+icen
	t_around(nk)=tout
	if(tout.lt.tmin) then
		index=1
		tmin=tout
		imin=icur+icen
	end if
	!write(*,*)icur+icen,' xrefl=',xrefl,' tout=',tout
134 continue
end do
if(index.eq.1) then
	icen=imin
	xrefl_min = xrefl0 + imin * dx_refl
	!write(*,*)' icen=',icen,' xrefl=',xrefl_min
	goto 135
end if

137 continue
dx_refl=dx_refl/2.
if(dx_refl.lt.dxmin_refl) goto 136
!write(*,*)' dx_refl=',dx_refl

do icase=1,2
	xrefl = xrefl_min + dx_refl*(-1)**icase
	yrefl = yyy_surf(xrefl)
	call streight_line(xzt,yzt,xrefl,yrefl,ips, tout1)
	call streight_line(xrefl,yrefl,xst,yst,ips, tout2)
	tout=tout1+tout2
	!write(*,*)icase,nodes,' xrefl=',xrefl,' tout=',tout
	if(tout.lt.tmin) then
		tmin=tout
		xrefl_min=xrefl
		!write(*,*)' icase=',icase,' xrefl=',xrefl_min,' tmin=',tmin
		exit
	end if
end do
goto 137

136 continue

xrefl = xrefl_min 
yrefl = yyy_surf(xrefl)

call streight_line(xzt,yzt,xrefl,yrefl,ips, tout1)

nrefl=nodes

ntmp=nodes
xtmp=xray
ytmp=yray

call streight_line(xrefl,yrefl,xst,yst,ips, tout2)

do i=2,nodes
	ntmp=ntmp+1
	xtmp(ntmp)=xray(i)
	ytmp(ntmp)=yray(i)
end do

nodes=ntmp
xray=xtmp
yray=ytmp

tout=tout1+tout2
xrefl=xray(nrefl)
yrefl=yray(nrefl)
!write(*,*)' 2 lines: xrefl=',xrefl,' yrefl=',yrefl,' tout=',tout
return
end