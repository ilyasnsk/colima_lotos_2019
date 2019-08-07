function vrefmod(z,ips)

common/refmod/nref,zref(600),vref(600,2)

vrefmod=-900.


if(nref.eq.0) then
    write(*,*)' the reference model is not defined!'
    pause
end if

if(z.le.zref(1))then
    vrefmod=vref(1,ips)
else if(z.ge.zref(nref)) then
    vrefmod=vref(nref,ips)
else
    do i=1,nref-1
        z1=zref(i)
        z2=zref(i+1)
        if((z-z1)*(z-z2).gt.0.)cycle
        v1=vref(i,ips)
        v2=vref(i+1,ips)
        vrefmod=v1+((v2-v1)/(z2-z1))*(z-z1)
        exit
    end do
end if


return
end 