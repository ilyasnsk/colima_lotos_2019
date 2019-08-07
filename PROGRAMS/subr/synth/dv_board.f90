function dv_board(xx,yy,zz,ips0)
common/brd_dv/xbr1(2),xbr2(2),dxbr1(2),dxbr2(2),&
			  ybr1(2),ybr2(2),dybr1(2),dybr2(2),&
			  zbr1(2),zbr2(2),dzbr1(2),dzbr2(2),amp(4)

dv_board=0
ips=ips0
if(ips0.gt.2) ips=ips0-2
if((xx-xbr1(ips))*(xx-xbr2(ips)).gt.0.) return
if((yy-ybr1(ips))*(yy-ybr2(ips)).gt.0.) return
if((zz-zbr1(ips))*(zz-zbr2(ips)).gt.0.) return


nxx=(xbr2(ips)-xbr1(ips))/(dxbr1(ips)+dxbr2(ips))
nyy=(ybr2(ips)-ybr1(ips))/(dybr1(ips)+dybr2(ips))
nzz=(zbr2(ips)-zbr1(ips))/(dzbr1(ips)+dzbr2(ips))
if(nzz.eq.0.) nzz=1

do ixx=1,1000
	x1=xbr1(ips)+(ixx-1)*(dxbr1(ips)+dxbr2(ips))
	x2=xbr1(ips)+ ixx   *(dxbr1(ips)+dxbr2(ips))
	if((xx-x1)*(xx-x2).le.0.) exit
end do

if(xx.ge.x1+dxbr1(ips)) return

do iyy=1,1000
	y1=ybr1(ips)+(iyy-1)*(dybr1(ips)+dybr2(ips))
	y2=ybr1(ips)+ iyy   *(dybr1(ips)+dybr2(ips))
	if((yy-y1)*(yy-y2).le.0.) exit
end do

if(yy.ge.y1+dybr1(ips)) return

do izz=1,1000
	z1=zbr1(ips)+(izz-1)*(dzbr1(ips)+dzbr2(ips))
	z2=zbr1(ips)+ izz   *(dzbr1(ips)+dzbr2(ips))
	if((zz-z1)*(zz-z2).le.0.) exit
end do

if(zz.ge.z1+dzbr1(ips)) return

sign=(-1.)**(ixx+iyy+izz)

dv_board=sign*amp(ips0)

return
end