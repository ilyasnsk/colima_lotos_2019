function velocity(x,y,z,ips)
common/flat_sf/key_flat1
vel_min=0.5

rz=6371

depth = sph_flat(key_flat1,x,y,z)

v0=vrefmod(depth,ips)

dv_apr= 0.01 * vert_anom(x,y,z,ips) * v0

dv = anom_3D_xyz_lin_v(x,y,z,ips)


velocity = v0 + dv + dv_apr

if(velocity .lt. vel_min) velocity = vel_min


return
end