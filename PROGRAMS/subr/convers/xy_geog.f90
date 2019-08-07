subroutine xy_geog(x,y,z,fi,tet,h)

common/pi/pi,per
common/povorot/fi0,tet0,orient,sinal,cosal
rz=6371.

xx = x*cosal+y*sinal
yy =-x*sinal+y*cosal
zz = rz-z

X1 =  xx
Y1 = -yy*SIN(tet0*PER) + zz*COS(tet0*PER)
Z1 =  yy*COS(tet0*PER) + zz*SIN(tet0*PER)

if(x1.eq.0.) then
	fi=fi0
else
	IF(X1.GT.0..AND.Y1.GT.0.)PA=0.
	IF(X1.LT.0.)PA=-PI
	IF(X1.GT.0..AND.Y1.LT.0.)PA=2.*PI
	ff=atan(y1/x1)/per
	FI=fi0-(ff+PA/PER)+90.
end if

if(abs(fi+360-fi0).lt.abs(fi-fi0))fi=fi+360
if(abs(fi-360-fi0).lt.abs(fi-fi0))fi=fi-360

r=sqrt(x1*x1+y1*y1+z1*z1)
tet=ASIN(Z1/r)/PER
h=rz-r


return
end