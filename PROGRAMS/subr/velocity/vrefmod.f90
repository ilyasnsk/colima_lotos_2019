function vrefmod(z,ips)
real zref(600)
common/refmod/nref,href(600),vref(600,2)

zref=href
if(nref.eq.0) then
	write(*,*)' the reference model is not defined!'
	pause
end if

!write(*,*)' z=',z,hmod(1),hmod(nref)
if(z.le.zref(1))then
	vrefmod=vref(ips,1)
else if(z.ge.zref(nref)) then
	z1=zref(nref-1)
	z2=zref(nref)
	v1=vref(ips,nref-1)
	v2=vref(ips,nref)
	vrefmod=v1+((v2-v1)/(z2-z1))*(z-z1)
else
	do i=1,nref-1
		z1=zref(i)
		z2=zref(i+1)
		v1=vref(ips,i)
		v2=vref(ips,i+1)
		if((z1-z)*(z2-z).gt.0.) cycle
		vrefmod=v1+((v2-v1)/(z2-z1))*(z-z1)
		exit
	end do
end if

return
end 