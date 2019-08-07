function dv_mod_2d(x,y,z,ips)
common/mod_2d/xmod1,nxmod,dxmod,ymod1,nymod,dymod, dv_mod(2,300,300,300)

dv_mod_2d=0
if(nxmod.eq.0) return

xmod2=xmod1+dxmod*nxmod
ymod2=ymod1+dymod*nymod
zmod2=zmod1+dzmod*nzmod

if((x-xmod1)*(x-xmod2).gt.0.) return
if((y-ymod1)*(y-ymod2).gt.0.) return
if((z-zmod1)*(z-zmod2).gt.0.) return

do ix=1,nxmod-1
	x1=xmod1+(ix-1)*dxmod
	x2=xmod1+(ix)*dxmod
	if((x-x1)*(x-x2).le.0.) exit
end do

do iy=1,nymod-1
	y1=ymod1+(iy-1)*dymod
	y2=ymod1+(iy)*dymod
	if((y-y1)*(y-y2).le.0.) exit
end do

do iz=1,nzmod-1
	z1=zmod1+(iz-1)*dzmod
	z2=zmod1+(iz)*dzmod
	if((z-z1)*(z-z2).le.0.) exit
end do

v111=dv_mod(ips,ix,iy,iz)
v211=dv_mod(ips,ix+1,iy,iz)
v121=dv_mod(ips,ix,iy+1,iz)
v221=dv_mod(ips,ix+1,iy+1,iz)

v112=dv_mod(ips,ix,iy,iz+1)
v212=dv_mod(ips,ix+1,iy,iz+1)
v122=dv_mod(ips,ix,iy+1,iz+1)
v222=dv_mod(ips,ix+1,iy+1,iz+1)

v11=v111+((v211-v111)/(x2-x1))*(x-x1)
v21=v121+((v221-v121)/(x2-x1))*(x-x1)
v12=v112+((v212-v112)/(x2-x1))*(x-x1)
v22=v122+((v222-v122)/(x2-x1))*(x-x1)

v1=v11+((v21-v11)/(y2-y1))*(y-y1)
v2=v12+((v22-v12)/(y2-y1))*(y-y1)

dv_mod_2d = v1+((v2-v1)/(z2-z1))*(z-z1)


return
end 