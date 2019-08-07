subroutine part_bending_z(ips,val_bend, tout)

common/ray_part/ nodes,xray(1000),yray(1000),zray(1000),sray(1000),kod_tmp(1000)

real xtmp(1000),ytmp(1000),ztmp(1000)
common/shift/ dxray(1000),dyray(1000),dzray(1000)

scur=0
do inode=2,nodes
	s1=sray(inode-1) 
	s2=sray(inode) 
	ds=s2-s1
	scur=scur+ds
end do

sss2=scur/2

sss=-sss2
dxray=0
dyray=0
dzray=0
scur=0
if(abs(val_bend).lt.0.000001) goto 432
do inode=2,nodes-1


	x1=xray(inode-1)
	y1=yray(inode-1)
	z1=zray(inode-1) 
	s1=sray(inode-1) 

	x2=xray(inode)
	y2=yray(inode)
	z2=zray(inode)
	s2=sray(inode)

	ds=s2-s1
	scur=scur+ds


	x3=xray(inode+1)
	y3=yray(inode+1)
	z3=zray(inode+1)

	dh=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1))
	dist=sqrt((z3-z1)*(z3-z1)+dh*dh)
	!write(*,*)' dh=',dh,' dist=',dist
!	write(*,*)' sss=',sss,' scur=',scur,' d=',d_value
	cosa=dh/dist

	if(dh.lt.0.0001) then
		cosb=1
		sinb=0
		sina=1
	else
		cosb=(x3-x1)/dh
		sinb=(y3-y1)/dh
		sina=(z3-z1)/dist
	end if

	sss=sss+ds

	d_value = val_bend * (1-(sss*sss)/(sss2*sss2))


	dxray(inode)= - d_value * sina * cosb
	dyray(inode)= - d_value * sina * sinb
	dzray(inode)= d_value * cosa
	!write(*,*)inode,dxray(inode),dyray(inode),dzray(inode)

end do
432 continue

xtmp = xray + dxray
ytmp = yray + dyray
ztmp = zray + dzray


!do i=1,nodes
!	write(*,*)xtmp(i),ytmp(i),ztmp(i)
!end do

ttt=0
sss=0
!write(*,*)' start'
do inode=1,nodes-1
	x1=xtmp(inode)
	x2=xtmp(inode+1)
	xm=(x1+x2)/2.

	y1=ytmp(inode)
	y2=ytmp(inode+1)
	ym=(y1+y2)/2.

	
	z1=ztmp(inode) 
	z2=ztmp(inode+1) 
	zm=(z1+z2)/2.

	!s1=sray(inode) 
	!s2=sray(inode+1) 
	!write(*,*)x1,x2
	!write(*,*)y1,y2
	!write(*,*)z1,z2

	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
	!ds2=s2-s1

	sss=sss+ds
	vvv=velocity (xm,ym,zm,ips)
	ttt=ttt+ds/vvv
	!write(*,*)' t=',ttt,' zm=',zm,' s=',sss,' v=',vvv
end do
!pause

tout=ttt

!write(*,*)nodes,' zmax=',zmax,' ttt=',ttt,' sss=',sss


return
end