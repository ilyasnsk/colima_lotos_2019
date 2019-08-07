function vert_anom(xx,yy,zzz,ips0)
common/dv_an/nan,anom(40,400,2),val(40,4),yan1(40),yan2(40),nnod(40)
common/profile/xa0(40),ya0(40),xb0(40),yb0(40)
common/center/fi0,tet0

vert_anom = 0

!write(*,*)' nan=',nan
if(nan.eq.0) return

ips=ips0; if(ips0.gt.2) ips=ips0-2

do ian=1,nan
	xa=xa0(ian)
	xb=xb0(ian)
	ya=ya0(ian)
	yb=yb0(ian)
	!write(*,*)' xa=',xa,' ya=',ya
	!write(*,*)' xb=',xb,' yb=',yb

	dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
	sinpov=(yb-ya)/dist
	cospov=(xb-xa)/dist

	xxx=(xx-xa)*cospov+(yy-ya)*sinpov
	yyy=-(xx-xa)*sinpov+(yy-ya)*cospov
	!write(*,*)xxx,yyy
	y1=yan1(ian)
	y2=yan2(ian)
	if((yyy-y1)*(yyy-y2).gt.0) cycle
	icrup=0

	do inod=2,nnod(ian)
		x1=anom(ian,inod-1,1)
		x2=anom(ian,inod,1)
		z1=anom(ian,inod-1,2)
		z2=anom(ian,inod,2)
        !write(*,*)x1,x2,xxx
		if(x1.eq.x2) cycle
		if ((xxx-x1)*(xxx-x2).gt.0.) cycle
		if ((xxx-x1).eq.0.) then
			in=inod
667			continue
			if(in.gt.2) then
				x0=anom(ian,in-2,1)
				if(x0.eq.x1) then
					in=in-1
					goto 667
				end if
			else
				x0=anom(ian,nnod(ian)-1,1)
			end if
			if ((x0-x1)*(x2-x1).gt.0.) cycle
		end if
		if ((xxx-x2).eq.0.) cycle
		zz=z1+((z2-z1)/(x2-x1))*(xxx-x1)
		!write(*,*)' xxx=',xxx,' zzz=',zzz,' zz=',zz
		if (zz.lt.zzz) cycle
		icrup=icrup+1
	end do
	icrup2=int(icrup/2)*2
!write(*,*)' ian=',ian,' icrup=',icrup,' icrup2=',icrup2
	if(icrup.eq.icrup2) cycle
	vert_anom=val(ian,ips0)
end do

return
end


