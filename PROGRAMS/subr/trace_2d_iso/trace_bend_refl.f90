subroutine trace_bend_refl(xzt,yzt,xst,yst,ips,	tout,xrefl,yrefl)

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



!************************************************************
!************************************************************
!************************************************************
!************************************************************



totdist=0
do i=1,nodes-1
	x1=xray(i)
	y1=yray(i)
	x2=xray(i+1)
	y2=yray(i+1)
	ds=sqrt((x2-x1)**2+(y2-y1)**2)
	totdist=totdist+ds
end do

!write(*,*)' totdist=',totdist,' nodes=',nodes,' tout=',tout

nsegm_max = int_best(totdist / ds_segm_min)

if(nsegm_max.eq.0) then
	return
end if

tmin=tout

do nsegm=1,nsegm_max

!open(22,file='nodes.dat')
	
	indexx=0
	dpart=1./nsegm
	val_cur0=bend_max0/nsegm
	xtmp=xray
	ytmp=yray

	do iseg=1,nsegm

		part1=(iseg-1)*dpart
		part2=iseg*dpart

		call part_ray(xzt,yzt,xst,yst,part1,part2)

		irfl_cur=0
		do inode=2,npart-1
			if(kod_part(inode).eq.nrefl)irfl_cur=inode
		end do
		!write(*,*)' iseg=',iseg,' irfl_cur=',irfl_cur



!write(22,*) xpart(1),ypart(1)		
		ttt=0
		do inode=1,npart-1
			x1=xpart(inode)
			x2=xpart(inode+1)
			xm=(x1+x2)/2.
			y1=ypart(inode)
			y2=ypart(inode+1)
			ym=(y1+y2)/2.
			ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
			vvv=velocity (xm,ym,ips)
			ttt=ttt+ds/vvv
			!write(*,*)inode,ypart(inode),dypart(inode),' v=',vvv
		end do
		tini=ttt
		!write(*,*)' tini=',tini,' npart=',npart

		dxpart=0
		dypart=0
		dxtmp=0
		dytmp=0
		valtot=0
		tmin=tini
		t_segm(iseg)=tmin
		val_cur=val_cur0
		

332		continue
		do icase=1,2
			ind=0
			val_bend = (-1)**icase * val_cur

331			continue
			val = val_bend 

			sss2=spart(npart)/2.
			sss=0

			if(irfl_cur.ne.0) then
				srefl=spart(irfl_cur)

				d_value = val * (1-abs(srefl-sss2)/sss2)
				!write(*,*)' srefl=',srefl,' sss2=',sss2
				!write(*,*)' val=',val,' d_value=',d_value
				xrefl_n=xpart(irfl_cur)+ d_value
				yrefl_n=yyy_surf(xrefl_n)
				yrefl=yyy_surf(xpart(irfl_cur))
				!write(*,*)' xrefl_n=',xrefl_n,' yrefl_n=',yrefl_n
				dx_refl=d_value
				dy_refl=yrefl_n-ypart(irfl_cur)
				!write(*,*)' dx_refl=',dx_refl,' dy_refl=',dy_refl

				dydx=dy_refl/dx_refl
				do inode=2,npart-1
					sss=spart(inode)
					d_value = val * (1-abs(sss-sss2)/sss2)
					dxpart(inode)=   d_value 
					dypart(inode)=   d_value * dydx
					!write(*,*)inode,dxtmp(inode),dytmp(inode)
				end do

			else
				do inode=2,npart-1
					
					if(kod_part(inode).eq.nrefl)irfl_cur=inode

					x1=xpart(inode-1)
					y1=ypart(inode-1)
						
					x2=xpart(inode)
					y2=ypart(inode)

					x3=xpart(inode+1)
					y3=ypart(inode+1)

					sss=spart(inode)

					dx=x3-x1
					dy=y3-y1
					ds2=sqrt(dx*dx+dy*dy)

					d_value = val * (1-abs(sss-sss2)/sss2)

					dxpart(inode)= - d_value * dy / ds2
					dypart(inode)=   d_value * dx / ds2
					!write(*,*)inode,d_value,dxtmp(inode),dytmp(inode)

				end do

			end if

			ttt=0
			sss=0
			do inode=1,npart-1
				x1=	xpart(inode)+	dxpart(inode)
				x2=	xpart(inode+1)+	dxpart(inode+1)
				xm=(x1+x2)/2.

				y1= ypart(inode) +	dypart(inode)
				y2=	ypart(inode+1)+	dypart(inode+1)
				ym=(y1+y2)/2.

				ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

				sss=sss+ds
				vvv=velocity (xm,ym,ips)
				ttt=ttt+ds/vvv
				!write(*,*)inode,' x=',x1,x2,' y=',y1,y2
				!write(*,*)inode,ypart(inode),dypart(inode),' v=',vvv
			end do
			tout=ttt

			!write(*,*)' val_bend=',val_bend,' tout=',tout

			if(tout.lt.tmin) then
				ind=ind+1
				tmin=tout
				valtot=valtot+val_bend
				do ipp=1,npart
					kod=kod_part(ipp)
					if(kod.eq.0) cycle
					xray(kod)=xray(kod)+dxpart(ipp)
					yray(kod)=yray(kod)+dypart(ipp)
					xpart(ipp)=xpart(ipp)+dxpart(ipp)
					ypart(ipp)=ypart(ipp)+dypart(ipp)
				end do
				!write(*,*)' valtot=',valtot,' val_bend=',val_bend,' tmin=',tmin
				call remeshing(nrefl)
				goto 331
			else
				if(ind.ne.0) exit
			end if

		end do

		val_cur = val_cur / 2.
		if(abs(val_cur).gt.bend_min0) goto 332

		t_segm(iseg)=tmin
		xrefl=xray(nrefl)
		yrefl=yray(nrefl)
		yrefl2=yyy_surf(xrefl)
		!write(*,*)' xrefl=',xrefl,' yrefl=',yrefl,yrefl2
		!write(*,*)' valtot=',valtot,' tmin=',tmin

	end do


333	continue
	ttot=0
	do iseg=1,nsegm
		!write(*,*)' indexx(ip)=',indexx(ip)
		ttot=ttot+t_segm(iseg)
	end do


	call remeshing(nrefl)

	!open(11,file='rays_refl.dat')
	!do i=1,nodes
		!write(11,*)xray(i),yray(i)
	!end do
	!close(11)
	

	!write(*,*)' val_cur0=',val_cur0,' ttot=',ttot

	!open(21,file='ray_iter.bln')
	!write(21,*)nodes
	!do i=1,nodes
		!write(21,*)xray(i),yray(i)
	!end do
	!close(21)
	!pause


end do
tout=ttot
return
end