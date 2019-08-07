function anom_3D_xyz_lin_v(xxx,yyy,zzz,ips)

real fx3(3),fy3(3),fz3(3),xabc(3),yabc(3),zabc(3)
integer ixabc(3),iyabc(3),izabc(3)

common/mod_3D/dv_3D(2,300,300,150),xx1,nxx,dxx,yy1,nyy,dyy,zz1,nzz,dzz
common/zlimits/zmin,zmax


anom_3D_xyz_lin_v=0
if(nxx.eq.0) return

!write(*,*)' xxx=',xxx,' yyy=',yyy,' zzz=',zzz

xx2=xx1+(nxx-1)*dxx
yy2=yy1+(nyy-1)*dyy
zz2=zz1+(nzz-1)*dzz

if((xxx-xx1)*(xxx-xx2).ge.0) return
if((yyy-yy1)*(yyy-yy2).ge.0) return
if((zzz-zz1)*(zzz-zz2).ge.0) return

do ix=1,nxx-1
	x1=xx1+(ix-1)*dxx
	x2=xx1+ix*dxx
	if((xxx-x1)*(xxx-x2).le.0) goto 10
end do
write(*,*)' problem in anom_3D_xyz:'
write(*,*)' xx1=',xx1,' xx2=',xx2,' xxx=',xxx
pause
10 continue
!write(*,*)' ix=',ix

do iy=1,nyy-1
	y1=yy1+(iy-1)*dyy
	y2=yy1+iy*dyy
	if((yyy-y1)*(yyy-y2).le.0) goto 11
end do
write(*,*)' problem in anom_3D_xyz:'
write(*,*)' yy1=',yy1,' yy2=',yy2,' yyy=',yyy
pause
11 continue
!write(*,*)' iy=',iy

do iz=1,nzz-1
	z1=zz1+(iz-1)*dzz
	z2=zz1+iz*dzz
	if((zzz-z1)*(zzz-z2).le.0) goto 12
end do
write(*,*)' problem in anom_3D_xyz:'
write(*,*)' zzz=',zzz
pause
12 continue
!write(*,*)' iz=',iz

v111=dv_3D(ips,ix,iy,iz)
v211=dv_3D(ips,ix+1,iy,iz)

v121=dv_3D(ips,ix,iy+1,iz)
v221=dv_3D(ips,ix+1,iy+1,iz)

!write(*,*)' ix=',ix,' iy=',iy,' iz=',iz

!write(*,*)v111,v211,v121,v221

v112=dv_3D(ips,ix,iy,iz+1)
v212=dv_3D(ips,ix+1,iy,iz+1)
v122=dv_3D(ips,ix,iy+1,iz+1)
v222=dv_3D(ips,ix+1,iy+1,iz+1)
!write(*,*)v112,v212,v122,v222

v11=v111+((v111-v211)/(x1-x2))*(xxx-x1)
v21=v121+((v121-v221)/(x1-x2))*(xxx-x1)
v12=v112+((v112-v212)/(x1-x2))*(xxx-x1)
v22=v122+((v122-v222)/(x1-x2))*(xxx-x1)

v1=v11+((v11-v12)/(z1-z2))*(zzz-z1)
v2=v21+((v21-v22)/(z1-z2))*(zzz-z1)

!write(*,*)' v1=',v1,' v2=',v2

anom_3D_xyz_lin_v=v1+((v1-v2)/(y1-y2))*(yyy-y1)

return
end