subroutine geog_xy(fi,tet,h, x,y,z)
common/pi/pi,per
common/povorot/fi0,tet0,orient,sinal,cosal
rz=6371.

x1 = (rz-h)*SIN((fi-fi0)*per)*COS(tet*per) 
y1 = (rz-h)*COS((fi-fi0)*per)*COS(tet*per) 
z1 = (rz-h)*SIN(tet*per) 

xx = x1
yy = -y1*SIN(tet0*per)+z1*COS(tet0*per) 
zz = rz-(y1*COS(tet0*per)+z1*SIN(tet0*per)) 

x  = xx*cosal-yy*sinal
y  = xx*sinal+yy*cosal
z  = zz


return
end