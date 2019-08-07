function epic_dist(fi1,tet1, fi2,tet2)
common/pi/pi,per

rz=6371.

cosf1=cos(fi1*per)
cost1=cos(tet1*per)
sinf1=sin(fi1*per)
sint1=sin(tet1*per)

cosf2=cos(fi2*per)
cost2=cos(tet2*per)
sinf2=sin(fi2*per)
sint2=sin(tet2*per)

x1=cosf1*cost1
y1=sinf1*cost1
z1=sint1

x2=cosf2*cost2
y2=sinf2*cost2
z2=sint2

dist = rz*sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
if(dist.lt.50) then
	epic_dist = (dist/rz)/per
else
	cosa = x1*x2 + y1*y2 + z1*z2
	epic_dist = ACOS(cosa)/PER
end if

RETURN
END
