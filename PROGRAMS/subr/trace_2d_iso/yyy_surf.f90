function yyy_surf(xxx)
common/surface/nsurf,xsurf(100),ysurf(100)

if(xxx.lt.xsurf(1)) then
	isrf=1
	x1=xsurf(1)
	x2=xsurf(2)
else if(xxx.gt.xsurf(nsurf)) then
	isrf=nsurf-1
	x1=xsurf(nsurf-1)
	x2=xsurf(nsurf)
else
	do isrf=1,nsurf-1
		x1=xsurf(isrf)
		x2=xsurf(isrf+1)
		if((x1-xxx)*(x2-xxx).le.0) exit
	end do
end if

y1=ysurf(isrf)
y2=ysurf(isrf+1)

yyy_surf = y1+((y2-y1)/(x2-x1))*(xxx-x1)




return
end