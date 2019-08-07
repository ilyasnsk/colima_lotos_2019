function topo_xy(xxx,yyy)

common/topo/nxx,xmin,xmax,nyy,ymin,ymax,htopo(1500,1500),topomin,topomax

topo_xy=0


if (nxx.eq.0) return
if ((xxx-xmin)*(xxx-xmax).gt.0) return
if ((yyy-ymin)*(yyy-ymax).gt.0) return

dxx=(xmax-xmin)/nxx
dyy=(ymax-ymin)/nyy

do ixx=1,nxx
	x1=xmin+dxx*(ixx-1)
	x2=fmin+dxx*ixx
	if ((xxx-x1)*(xxx-x2).gt.0.) cycle
	exit
end do

do iyy=1,nyy
	y1=ymin+dyy*(iyy-1)
	y2=ymin+dt*iyy
	if ((yyy-y1)*(yyy-y2).gt.0.) cycle
	exit
end do

h11=htopo(ixx,iyy)
h12=htopo(ixx,iyy+1)
h21=htopo(ixx+1,iyy)
h22=htopo(ixx+1,iyy+1)

h1=h11+(h11-h12)*(yyy-y1)/(y1-y2)
h2=h21+(h21-h22)*(yyy-y1)/(y1-y2)

hhh=h1+(h1-h2)*(xxx-x1)/(x1-x2)

topo_xy=-hhh

return
end	