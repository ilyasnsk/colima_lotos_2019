function GAP_angle(fzt,tzt)

real azkr(5000),aztmp(5000)

common/rays_1ztr/ nkrat,istkr(5000),ipskr(5000)
common/stations/fstat(500),tstat(500),zstat(500)
common/pi/pi,per





do ikr=1,nkrat

	ist=istkr(ikr)
	fst=fstat(ist)
	tst=tstat(ist)

	dist=epic_dist(fzt,tzt, fst,tst)

	call SFDEC(fst,tst,0.,X,Y,Z,fzt,tzt)
	if(abs(x).gt.0.1) then
		az = atan(y/x)/per
		if(x.lt.0.) az=180+az
	else
		az = 180.
		if(y.lt.0.) az=-180
	end if
	azkr(ikr)=az
	!write(*,*)' az=',az
end do

aztmp=-999
nn=0
do ikr=1,nkrat
	azmin=900
	do ik2=1,nkrat
		if(azkr(ik2).gt.azmin) cycle
		azmin=azkr(ik2)
		imin=ik2
	end do
	nn=nn+1
	aztmp(nn)=azmin
	azkr(imin)=999
end do

freemax=0.
do i=1,nn
	if(i.eq.1) then
		a1=aztmp(nn)-360
	else
		a1=aztmp(i-1)
	end if
	a2=aztmp(i)
	da=a2-a1
!		write(*,*)a1,a2,da
	if(da.gt.freemax)freemax=da
end do

GAP_angle=freemax
!write(*,*)' freemax=',freemax,' nkrat=',nkrat


return
end