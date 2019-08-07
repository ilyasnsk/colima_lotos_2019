subroutine streight_line(xzt,yzt,zzt, xst,yst,zst, ips, tout)
common/ray/ nodes,xray(1000),yray(1000),zray(1000)
common/ray_param/ds_ini,ds_part_min,bend_min0,bend_max0


dist = sqrt((xst-xzt)*(xst-xzt)+(yst-yzt)*(yst-yzt)+(zst-zzt)*(zst-zzt))


nodes = int_best(dist/ds_ini)+1
ds = dist / (nodes-1)
!write(*,*)' dist=',dist,' nodes=',nodes,' ds=',ds

if(nodes.eq.1) then
	nodes=2
	ds=dist
end if

!write(*,*)' ztr:',xzt,yzt
do inode = 1,nodes
	sss = ds * (inode-1)
	xray(inode) = xzt + sss * (xst-xzt)/dist
	yray(inode) = yzt + sss * (yst-yzt)/dist
	zray(inode) = zzt + sss * (zst-zzt)/dist
	!write(*,*) sss,xray(inode),yray(inode)
end do
!write(*,*)sss,' sta:',xst,yst

ddd=0
ttt=0
sss=0
do inode=1,nodes-1
	x1=xray(inode)
	x2=xray(inode+1)
	xm=(x1+x2)/2.

	y1=yray(inode)
	y2=yray(inode+1)
	ym=(y1+y2)/2.

	z1=zray(inode)
	z2=zray(inode+1)
	zm=(z1+z2)/2.

	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
	!vvv=vrefmod (ym)
	vvv=velocity (xm,ym,zm,ips)
    !write(*,*)xm,ym,zm,vvv
	!write(*,*)' x 1 2=',x1,x2,' y 1 2=',y1,y2
	ttt=ttt+ds/vvv
	sss=sss+ds
end do
!write(*,*)' ttt=',ttt,' sss=',sss
!pause

tout=ttt


return
end