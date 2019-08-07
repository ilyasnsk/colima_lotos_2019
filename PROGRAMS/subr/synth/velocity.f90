function velocity(x,y,z,ips)
common/nanom/n_anomaly
common/center/fi0,tet0


v0=vrefmod(z,ips)

!call decsf(x,y,0.,fi0,tet0,fi,tet,h)

dv=anomaly(x,y,z,ips)

velocity = v0*(1 + 0.01 * dv)

!write(*,'(5f8.3)')x,y,z,v0,dv



!v0=vref_smooth(z,ips)
!dv = anom_3D_xyz(x,y,z,ips)


!velocity = v0 + dv


return
end