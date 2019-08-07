function alin(xcur, n,x0,dx,fun)
real fun(n)
do i=1,n-1
	x1=x0+dx*(i-0.5)
	if(i.eq.1)x1=x0
	x2=x0+dx*(i+0.5)
	if((x1-xcur)*(x2-xcur).gt.0.)cycle
	f1=fun(i)
	f2=fun(i+1)
	alin=f1+((f2-f1)/(x2-x1))*(xcur-x1)
	return
end do
alin=0.
return
end


function alin2d(x,y, x0,nx,dx,y0,ny,dy,fun)
real fun(nx,ny)
do ix=1,nx-1
	x1=x0+dx*(ix-0.5)
	if(ix.eq.1)x1=x0
	x2=x0+dx*(ix+0.5)
	if(ix.eq.nx-1)x2=x0+dx*(ix+1)
	if((x1-x)*(x2-x).gt.0.)cycle
	do iy=1,ny-1
		y1=y0+dy*(iy-0.5)
		if(iy.eq.1)y1=y0
		y2=y0+dy*(iy+0.5)
		if(iy.eq.ny-1)y2=y0+dy*(iy+1)
		if((y1-y)*(y2-y).gt.0.)cycle
		f11=fun(ix,iy)
		f12=fun(ix,iy+1)
		f21=fun(ix+1,iy)
		f22=fun(ix+1,iy+1)
		f1=f11+((f12-f11)/(y2-y1))*(y-y1)
		f2=f21+((f22-f21)/(y2-y1))*(y-y1)
		alin2d=f1+((f2-f1)/(x2-x1))*(x-x1)
		return
	end do
end do
alin=0.
return
end
