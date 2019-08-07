function vrefmod(z,ips)

common/refmod/nrefmod,zref(600),vref(600,2)

vrefmod=-900.


if(nrefmod.eq.0) then
    write(*,*)' z=',z,' ips=',ips
    write(*,*)' the reference model is not defined!'
    pause
end if

if(z.le.zref(1))then
    vrefmod=vref(1,ips)
else if(z.ge.zref(nrefmod)) then
    vrefmod=vref(nrefmod,ips)
else
    do i=1,nrefmod-1
        z1=zref(i)
        z2=zref(i+1)
        if((z-z1)*(z-z2).gt.0.)cycle
        v1=vref(i,ips)
        v2=vref(i+1,ips)
        !write(*,*)' i=',i,' v1=',v1,' v2=',v2
        vrefmod=v1+((v2-v1)/(z2-z1))*(z-z1)
        exit
    end do
end if


return
end 