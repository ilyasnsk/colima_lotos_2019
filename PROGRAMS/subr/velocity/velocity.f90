function velocity(x,y,z,ips)

v0=vrefmod(z,ips)



dv_prev=dv_mod_2d(x,y,z,ips)

velocity = v0 + dv_prev

!write(*,*)' velocity=',v0,dv_prev
if(velocity.lt.0.00001) stop

return
end