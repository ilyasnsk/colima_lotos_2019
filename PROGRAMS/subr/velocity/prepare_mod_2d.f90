subroutine prepare_mod_2d(md,iter)
character*8 md
character*1 it,ps
common/mod_2d/xmod1,nxmod,dxmod,ymod1,nymod,dymod, dv_mod(2,300,300)

nxmod=0
nymod=0

if(iter.eq.0) return

write(it,'(i1)')iter

do ips=1,2
	write(ps,'(i1)')ips
	open(1,file='../../DATA/'//md//'/data/dv_prev'//it//ps//'.dat',form='binary')
	read(1)xmod1,xmod2,dxmod,nxmod
	read(1)ymod1,ymod2,dymod,nymod

	write(*,*)' ips=',ips,' nxmod=',nxmod,' nymod=',nymod
	read(1)((dv_mod(ips,ix,iy),ix=1,nxmod),iy=1,nymod)
	close(1)
end do

return
end