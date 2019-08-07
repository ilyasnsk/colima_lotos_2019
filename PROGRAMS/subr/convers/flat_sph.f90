function flat_sph(key_flat1,x,y,depth)
rz=6371

if(key_flat1.eq.2) then
    zz = rz - sqrt( (rz-depth)*(rz-depth) - (x*x + y*y) )
else
    zz = depth
end if 

flat_sph = zz

return
end