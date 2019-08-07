function our_noise(nray,xmid,ips)

USE DFPORT

common/histo/istep,xconv(2000),xdir1,xdir2,dxdir
common/noise_1/a_rand(2),n_perc_out,a_out,a_stat

our_noise=0

nfreq_out=10000000
if(n_perc_out.ne.0)nfreq_out=100/n_perc_out


xmid1=0
x = RAND()
!write(*,*)' x=',x
do ii=1,istep-1
	x1=xdir1+dxdir*(ii-1)
	x2=xdir1+dxdir*ii
	if(x.lt.x1.or.x.ge.x2) cycle
	y1=xconv(ii)
	y2=xconv(ii+1)
	y=y1+((y2-y1)/(x2-x1))*(x-x1)
	!write(*,*)' y=',y
	exit
end do

a_noise=a_rand(ips)
if(mod(nray,nfreq_out).eq.0) a_noise = a_rand(ips) * a_out

our_noise = y * a_noise

return
end

