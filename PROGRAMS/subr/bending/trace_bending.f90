subroutine trace_bending(xzt,yzt,zzt,xst,yst,zst,ips,	tout)
common/ray/ nodes,xray(1000),yray(1000),zray(1000)
common/ray_min/ nrmin,xrmin(1000),yrmin(1000),zrmin(1000)
common/ray_part/ nrpart,xpart(1000),ypart(1000),zpart(1000),spart(1000),kod_part(1000)
common/shift/ dxray(1000),dyray(1000),dzray(1000)
common/ray_param/ds_ini,ds_part_min,val_bend_min,bend_max0


real xrmin2(1000),yrmin2(1000),zrmin2(1000)

real t_part(100)
integer indexx(100)

!ds_part_min = 5
!val_bend_min = 0.02
!bend_max0 = 10


call straight_line(xzt,yzt,zzt,xst,yst,zst,ips, tout)
totdist=sqrt((xzt-xst)*(xzt-xst)+(yzt-yst)*(yzt-yst)+(zzt-zst)*(zzt-zst))

!write(*,*)' straight line, tout=',tout

npart_max = int_best(totdist / ds_part_min)
val_bend0 = val_bend_min * npart_max
!write(*,*)' npart_max=',npart_max,' val_bend0=',val_bend0

if(npart_max.eq.0) then
	return
end if

tmin=tout
nrmin=nodes
xrmin=xray
yrmin=yray
zrmin=zray

!goto 12

sumbend=0

bend=bend_max0/10

11 continue

	sumbend = sumbend + bend

	call part_ray(xzt,yzt,zzt,xst,yst,zst,ips,0.,1.)


	!write(*,*)' nodes=',nodes,' nrpart=',nrpart
	!call part_bending_z(ips,0., tini)
	!write(*,*)' tini=',tini

	call part_bending_z(ips,bend, tout)

	nodes=nrpart
	xray=xpart+dxray
	yray=ypart+dyray
	zray=zpart+dzray

	!call remeshing()
	zmax=-999.
	do i=1,nodes
		!write(*,*)' zray(i)=',zray(i),' zmax=',zmax
		if(zray(i).gt.zmax) zmax=zray(i)
	end do
	!write(*,*)nodes,' sum=',sumbend,' zmax=',zmax,' t=',tout

	if(tout.lt.tmin) then
		tmin=tout
		nrmin=nodes
		xrmin=xray
		yrmin=yray
		zrmin=zray
		goto 11
	end if

!if(sumbend.lt.bend_max0) goto 11

12 continue

do nparts=1,npart_max

243	continue
	indexx=0
	dpart=1./nparts
	nodes=nrmin
	!do i=1,nodes
		!write(*,*)zray(i),zrmin(i),zray(i)-zrmin(i)
	!end do
	xray=xrmin
	yray=yrmin
	zray=zrmin



	do ip=1,nparts

		!if(indexx(ip).eq.0) cycle

		part1=(ip-1)*dpart
		part2=ip*dpart


		!call remeshing()

		call part_ray(xzt,yzt,zzt,xst,yst,zst,ips,part1,part2)

		!do inode=1,nrpart
		!	write(*,*)kod_part(inode),spart(inode),xpart(inode),ypart(inode),zpart(inode)
		!end do

!write(*,*)' nodes=',nodes,' nrpart=',nrpart
		call part_bending_hor(ips,0., tini)

!write(*,*)' part=',ip,' tini=',tini,' s=',spart(nrpart)

		tmin=tini
		t_part(ip)=tmin
		!indexx(ip)=0
		!cycle

! Vertical bending:
		do icase=1,2
			val_bend=-(-1)**icase*(val_bend0/nparts)
			ind=1

331			continue
			val = val_bend * ind
			call part_bending_z(ips,val, tout)
			!write(*,*)' ind=',ind,' val_bend=',val_bend
			!write(*,*)' ver shift: ip=',ip,' val=',val,' dt=',tout-tmin

			if(tout.lt.tmin) then
				tmin=tout
				!xrmin=xray
				!yrmin=yray
				!zrmin=zray
				do ipp=1,nrpart
					kod=kod_part(ipp)
					!write(*,*)ipp,kod,dxray(ipp),dyray(ipp),dzray(ipp)
					if(kod.eq.0) cycle
					xrmin(kod)=xray(kod)+dxray(ipp)
					yrmin(kod)=yray(kod)+dyray(ipp)
					zrmin(kod)=zray(kod)+dzray(ipp)
					!if(nparts.eq.9)write(*,*)kod,dxray(ipp),dyray(ipp),dyray(ipp)
				end do
				ind=ind+1
				goto 331	! Try a larger shift
			else
				if(ind.gt.1) then
					nodes=nrmin
					xray=xrmin
					yray=yrmin
					zray=zrmin
					indexx(ip)=1
					t_part(ip)=tmin
					call part_ray(xzt,yzt,zzt,xst,yst,zst,ips,part1,part2)
					exit
				end if 
			end if
		end do
		!write(*,*)' tmin=',tmin

! Horizontal bending:
		do icase=1,2
			val_bend=(-1)**icase*(val_bend0/nparts)
			ind=1


334			continue
			val = val_bend * ind
!if(ip.eq.15.and.nparts.eq.16) then
	!write(*,*)' val=',val,' tmin=',tmin
	!call part_bending_hor_7(ips,val, tout)
	!pause
!else
			call part_bending_hor(ips,val, tout)
!end if
			if(tout.lt.tmin) then
				!write(*,*)ip,' val=',val,' dt=',tout-tmin,' t=',tout
				tmin=tout
				!xrmin=xray
				!yrmin=yray
				!zrmin=zray
				do ipp=1,nrpart
					kod=kod_part(ipp)
					if(kod.eq.0) cycle
					xrmin(kod)=xray(kod)+dxray(ipp)
					yrmin(kod)=yray(kod)+dyray(ipp)
					zrmin(kod)=zray(kod)+dzray(ipp)
!if(ip.eq.15.and.nparts.eq.16)write(*,*)kod,' yray:',yray(kod),dyray(ipp),' ymin:',yrmin(kod)
				end do
				ind=ind+1
				goto 334
			else
				if(ind.gt.1) then
					nodes=nrmin
					xray=xrmin
					yray=yrmin
					zray=zrmin
					indexx(ip)=1
					t_part(ip)=tmin
					call part_ray(xzt,yzt,zzt,xst,yst,zst,ips,part1,part2)
					exit
				end if 
			end if
		end do		! 2 attempts 

	end do
	ind_tot=0
	ttot=0
	do ip=1,nparts
		!write(*,*)' indexx(ip)=',indexx(ip)
		ttot=ttot+t_part(ip)
		if(indexx(ip).eq.1) ind_tot=1
	end do
	!write(*,*)' nparts=',nparts,' ttot=',ttot
	!if(ind_tot.eq.1) goto 243

end do
tout=ttot
!write(*,*)' tout=',tout

return

end
