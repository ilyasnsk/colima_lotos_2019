function sph_flat(key_flat1,x,y,z)
rz=6371

if(key_flat1.eq.2) then
    depth = rz - sqrt( (rz-z)*(rz-z) + (x*x + y*y) )
else
    depth = z 
end if 

sph_flat = depth

return
end