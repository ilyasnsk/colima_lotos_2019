function relief_surf(fi,tet)

common/topo/hsurf,nfi,fmin,fmax,ntet,tmin,tmax,htopo(1500,1500)

hhh=0.

!write(*,*)' nfi=',nfi,' fi=',fi,' fmin=',fmin,' fmax=',fmax

if (nfi.eq.0) then
	h_lim=0
	return
end if

if ((fi-fmin)*(fi-fmax).gt.0) then
	h_lim=0
	return
end if

if ((tet-tmin)*(tet-tmax).gt.0) then
	h_lim=0
	return
end if

df=(fmax-fmin)/nfi
dt=(tmax-tmin)/ntet

do ifi=1,nfi
	f1=fmin+df*(ifi-1)
	f2=fmin+df*ifi
	if ((fi-f1)*(fi-f2).gt.0.) cycle
	exit
end do

do itet=1,ntet
	t1=tmin+dt*(itet-1)
	t2=tmin+dt*itet
	if ((tet-t1)*(tet-t2).gt.0.) cycle
	exit
end do

h11=htopo(ifi,itet)
h12=htopo(ifi,itet+1)
h21=htopo(ifi+1,itet)
h22=htopo(ifi+1,itet+1)


h1=h11+(h11-h12)*(tet-t1)/(t1-t2)
h2=h21+(h21-h22)*(tet-t1)/(t1-t2)

hhh=h1+(h1-h2)*(fi-f1)/(f1-f2)

relief_surf=-hhh

return
end	