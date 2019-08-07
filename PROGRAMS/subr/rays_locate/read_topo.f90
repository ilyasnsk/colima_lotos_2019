subroutine read_topo(ar)
character*8 md,ar
common/topo/hsurf,nfi,fmin,fmax,ntet,tmin,tmax,htopo(1500,1500)

nfi=0
hsurf=0.
open(1,file='../../../DATA/'//ar//'/inidata/topo.grd',status='old',err=1)

read(1,*)
read(1,*)nfi,ntet
read(1,*)fmin,fmax
read(1,*)tmin,tmax
read(1,*)hmin,hmax
do itet=1,ntet
	read(1,*)(htopo(ifi,itet),ifi=1,nfi)
end do
close(1)

do ifi=1,nfi
	do itet=1,ntet
		htopo(ifi,itet)=htopo(ifi,itet)/1000.
	end do
end do
!write(*,*)' nar=',nar

1 continue
return
end