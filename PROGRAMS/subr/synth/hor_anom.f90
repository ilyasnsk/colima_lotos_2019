function hor_anom(fi,tet,zzz,ips0)

common/dv_hor_an/nan,anom(40,400,2),val(40,4),zan1(40),zan2(40),nnod(40)


!write(*,*)' nan=',nan

ips=ips0
if(ips0.gt.2) ips=ips0-2
xxx=fi
yyy=tet


hor_anom = 0

if(nan.eq.0) return



do ian=1,nan
	z1=zan1(ian)
	z2=zan2(ian)
	!write(*,*)'nnod=',nnod(ian),' z1=',z1,' z2=',z2
	if((zzz-z1)*(zzz-z2).gt.0) cycle
	icrup=0

	do inod=2,nnod(ian)
!write(*,*)' fi=',anom(ian,inod,1),' tet=',anom(ian,inod,2)
		x1=anom(ian,inod-1,1)
		x2=anom(ian,inod,1)
		y1=anom(ian,inod-1,2)
		y2=anom(ian,inod,2)
!write(*,*)' xxx=',xxx,' x1=',x1,' x2=',x2
!write(*,*)' yyy=',yyy,' y1=',y1,' y2=',y2
!pause
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
		yy=y1+((y2-y1)/(x2-x1))*(xxx-x1)
		!write(*,*)' x1=',x1,' x2=',x2,' yy=',yy
		if (yy.lt.yyy) cycle
		icrup=icrup+1
	end do
	icrup2=int(icrup/2)*2
!write(*,*)' ian=',ian,' icrup=',icrup,' icrup2=',icrup2
	if(icrup.eq.icrup2) cycle
	hor_anom=val(ian,ips0)
end do

return
end


