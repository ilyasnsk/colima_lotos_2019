function our_noise(ips)

USE DFPORT

common/histo/istep,xconv(2000),xdir1,xdir2,dxdir
common/noise_1/a_rand(2),xmid(2)

our_noise=0

x = RAND()
do ii=1,istep-1
	x1=xdir1+dxdir*(ii-1)
	x2=xdir1+dxdir*ii
	if(x.lt.x1.or.x.ge.x2) cycle
	y1=xconv(ii)
	y2=xconv(ii+1)
	y=y1+((y2-y1)/(x2-x1))*(x-x1)
	!write(*,*)' y=',y
	goto 16
end do

write(*,*)' x=',x,' istep=',istep,' dxdir=',dxdir
write(*,*)' xconv(1)=',xconv(1),' xconv(istep)=',xconv(istep)
write(*,*)' xdir1=',xdir1,' xdir2=',xdir2
stop

16 continue

our_noise = xmid(ips) + y * a_rand(ips)

return
end

