function vert_brd(xx,yy,zzz,ips0)
common/brd_vert/nprof,fia(10),teta(10),fib(10),tetb(10),ybr1(10),ybr2(10),amp(10,4),&
				xbr1(10,2),xbr2(10,2),dxbr1(10,2),dxbr2(10,2),&
				zbr1(10,2),zbr2(10,2),dzbr1(10,2),dzbr2(10,2)


common/center/fi0,tet0

ips=ips0
if(ips0.gt.2) ips=ips0-2
vert_brd = 0

do ipr=1,nprof

	call SFDEC(fia(ipr),teta(ipr),0.,xa,ya,Z,fi0,tet0)
	call SFDEC(fib(ipr),tetb(ipr),0.,xb,yb,Z,fi0,tet0)

	dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
	sinpov=(yb-ya)/dist
	cospov=(xb-xa)/dist

	xxx=(xx-xa)*cospov+(yy-ya)*sinpov
	yyy=-(xx-xa)*sinpov+(yy-ya)*cospov

	!write(*,*)' xxx=',xxx,' yyy=',yyy,' zzz=',zzz


	if((yyy-ybr1(ipr))*(yyy-ybr2(ipr)).gt.0) cycle

	if((xxx-xbr1(ipr,ips))*(xxx-xbr2(ipr,ips)).gt.0.) cycle
	if((zzz-zbr1(ipr,ips))*(zzz-zbr2(ipr,ips)).gt.0.) cycle

	nx=(xbr2(ipr,ips)-xbr1(ipr,ips))/(dxbr1(ipr,ips)+dxbr2(ipr,ips))
	nz=(zbr2(ipr,ips)-zbr1(ipr,ips))/(dzbr1(ipr,ips)+dzbr2(ipr,ips))
	if(nz.eq.0.) nz=1

	do ix=1,nx
		x1=xbr1(ipr,ips)+(ix-1)*(dxbr1(ipr,ips)+dxbr2(ipr,ips))
		x2=xbr1(ipr,ips)+ ix   *(dxbr1(ipr,ips)+dxbr2(ipr,ips))
		if((xxx-x1)*(xxx-x2).le.0.) exit
	end do
	if(xxx.ge.x1+dxbr1(ipr,ips)) cycle

	do iz=1,nz
		z1=zbr1(ipr,ips)+(iz-1)*(dzbr1(ipr,ips)+dzbr2(ipr,ips))
		z2=zbr1(ipr,ips)+ iz   *(dzbr1(ipr,ips)+dzbr2(ipr,ips))
		if((zzz-z1)*(zzz-z2).le.0.) exit
	end do
	if(zzz.ge.z1+dzbr1(ipr,ips)) cycle

	sign=(-1.)**(ix+iz)

	vert_brd=sign*amp(ipr,ips0)

end do




return
end


