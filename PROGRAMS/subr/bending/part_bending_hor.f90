subroutine part_bending_hor(ips,val_bend, tout)

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


xa=xray(1)
xb=xray(nodes)
ya=yray(1)
yb=yray(nodes)
dh=sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb))
cosb=(xa-xb)/dh
sinb=(ya-yb)/dh

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

	sss=sss+ds

	d_value = val_bend * (1-(sss*sss)/(sss2*sss2))

	!write(*,*)' sss=',sss,' scur=',scur,' d=',d_value

	dxray(inode)= - d_value * sinb
	dyray(inode)=   d_value * cosb
	dzray(inode)= 0.

	!write(*,*)inode,dxray(inode),dyray(inode),dzray(inode)

end do
432 continue

xtmp = xray + dxray
ytmp = yray + dyray
ztmp = zray + dzray


!do i=1,nodes
!	write(*,*)xray(i),yray(i),zray(i)
!end do

ttt=0
sss=0
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


	s1=sray(inode) 
	s2=sray(inode+1) 

	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
	ds2=s2-s1
	!write(*,*)inode,ds,ds2

	sss=sss+ds
	vvv=velocity (xm,ym,zm,ips)
	ttt=ttt+ds/vvv
	!write(*,*)inode,xm,ym,zm,vvv
end do

tout=ttt

!write(*,*)nodes,' zmax=',zmax,' ttt=',ttt,' sss=',sss


return
end