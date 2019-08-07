subroutine trace_bend_2d(xzt,yzt,xst,yst,ips,	tout)

real xrmin(1000),yrmin(1000)
real t_segm(1000)
integer indexx(100)
real dxtmp(1000),dytmp(1000),xtmp(1000),ytmp(1000)

common/ray_param/ds_ini,ds_segm_min,bend_min0,bend_max0

common/ray/ nodes,xray(1000),yray(1000)
common/ray_part/ npart,xpart(1000),ypart(1000),spart(1000),kod_part(1000)
common/shift/ dxpart(1000),dypart(1000)

totdist=sqrt((xzt-xst)*(xzt-xst)+(yzt-yst)*(yzt-yst))
tout=0

if(totdist.lt.ds_segm_min) then
		x1=xzt
		x2=xst
		xm=(x1+x2)/2.
		y1=yzt
		y2=yst
		ym=(y1+y2)/2.

		vvv=velocity (xm,ym,ips)
		tout=ds/vvv
		return
end if
nsegm_max = int_best(totdist / ds_segm_min)

if(nsegm_max.eq.0) then
	return
end if


call streight_line(xzt,yzt,xst,yst,ips, tout)
!write(*,*)' streight line: tout=',tout
!write(*,*)ds_ini,ds_segm_min,bend_min0,bend_max0
!pause

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
		!write(*,*)' iseg=',iseg,' part1=',part1,' part2=',part2

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
			do inode=2,npart-1

				x1=xpart(inode-1)
				y1=ypart(inode-1)
				x2=xpart(inode)
				y2=ypart(inode)
				x3=xpart(inode+1)
				y3=ypart(inode+1)

				ds=sqrt((x2-x1)**2+(y2-y1)**2)
				sss=sss+ds

				dx=x3-x1
				dy=y3-y1
				ds2=sqrt(dx*dx+dy*dy)

				d_value = val * (1-abs(sss-sss2)/sss2)

				dxtmp(inode)= - d_value * dy / ds2
				dytmp(inode)=   d_value * dx / ds2
				!write(*,*)inode,d_value,dxpart(inode),dypart(inode)

			end do

			ttt=0
			sss=0
			do inode=1,npart-1
				x1=	xpart(inode)+	dxpart(inode)+		dxtmp(inode)
				x2=	xpart(inode+1)+	dxpart(inode+1)+	dxtmp(inode+1)
				xm=(x1+x2)/2.

				y1=	ypart(inode)+	dypart(inode)+		dytmp(inode)
				y2=	ypart(inode+1)+	dypart(inode+1)+	dytmp(inode+1)
				ym=(y1+y2)/2.

				ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

				sss=sss+ds
				vvv=velocity (xm,ym,ips)
				ttt=ttt+ds/vvv
				!write(*,*)inode,ypart(inode),dypart(inode),' v=',vvv
			end do

			tout=ttt
			!write(*,*)' val_bend=',val_bend,' tout=',tout

			if(tout.lt.tmin) then
				ind=ind+1
				tmin=tout
				dxpart=dxpart+dxtmp
				dypart=dypart+dytmp
				valtot=valtot+val_bend
				!write(*,*)' valtot=',valtot,' val_bend=',val_bend,' tmin=',tmin
				goto 331
			else
				if(ind.ne.0) exit
			end if

		end do

		val_cur = val_cur / 2.
		if(abs(val_cur).gt.bend_min0) goto 332


		do ipp=1,npart
			kod=kod_part(ipp)
			if(kod.eq.0) cycle
			xtmp(kod)=xtmp(kod)+dxpart(ipp)
			ytmp(kod)=ytmp(kod)+dypart(ipp)
			xpart(ipp)=xpart(ipp)+dxpart(ipp)
			ypart(ipp)=ypart(ipp)+dypart(ipp)
		end do
		t_segm(iseg)=tmin

	end do


333	continue
	ttot=0
	do iseg=1,nsegm
		!write(*,*)' indexx(ip)=',indexx(ip)
		ttot=ttot+t_segm(iseg)
	end do

!write(22,*)xray(nodes),yray(nodes)
!close(22)

	xray=xtmp
	yray=ytmp
	!write(*,*)' val_cur0=',val_cur0,' ttot=',ttot
	
	nrefl=0
	call remeshing(nrefl)


	!open(21,file='ray_iter.bln')
	!write(21,*)nodes
	!do i=1,nodes
		!write(21,*)xray(i),yray(i)
	!end do
	!close(21)
	!pause


end do
tout=ttot

ttt=0
sss=0
do inode=1,nodes-1
	x1=	xray(inode)
	x2=	xray(inode+1)
	xm=(x1+x2)/2.

	y1=	yray(inode)
	y2=	yray(inode+1)
	ym=(y1+y2)/2.

	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

	sss=sss+ds

	vvv=velocity (xm,ym,ips)
	!dv= dv_mod_2d(xm,ym)

!	write(*,*)xm,ym,' vvv=',vvv
	!call vel_mod_2d(xm,ym, dv000,dv060,dv120)
	!write(*,*)' dv000=',dv000,' dv060=',dv060,' dv120=',dv120


	ttt=ttt+ds/vvv
	!write(*,*)inode,ypart(inode),dypart(inode),' v=',vvv
end do

tout=ttt




return
end