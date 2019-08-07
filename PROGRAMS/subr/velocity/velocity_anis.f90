function velocity(x,y)

v0=vrefmod(y)

dv_prev=dv_mod_2d(x,y)

velocity = v0 + dv_prev


return
end