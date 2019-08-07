subroutine read_3D_mod_v(ar,md,iter)
character*8 ar,md
character*1 ps,rm,it 
common/mod_3D/dv_3D(2,300,300,150),xx1,nxx,dxx,yy1,nyy,dyy,zz1,nzz,dzz
common/zlimits/zmin,zmax
nxmax=300
nymax=300
nzmax=150

nxx=0
nyy=0
nzz=0

if(iter.lt.1) return

write(it,'(i1)')iter

dv_3D=0

do ips=1,2
	write(ps,'(i1)')ips
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/dv_v'//ps//it//'.dat',form='binary')

	read(1,end=13)xx10,nxx0,dxx0
	!write(*,*)xx10,nxx0,dxx0
	if(nxx0.eq.0) cycle
	if(nxx0.gt.nxmax) then
		write(*,*)' nxx > nxmax!'
		write(*,*)' Value of the nxx in common file "mod_3D" should be increased'
		pause
	end if
	xx1=xx10
	nxx=nxx0
	dxx=dxx0

	read(1)yy1,nyy,dyy
	if(nyy.gt.nymax) then
		write(*,*)' nyy > nymax!'
		write(*,*)' Value of the nyy in common file "mod_3D" should be increased'
		pause
	end if


	read(1)zz1,nzz,dzz
	if(nzz.gt.nzmax) then
		write(*,*)' nzz > nzmax!'
		write(*,*)' Value of the nzz in common file "mod_3D" should be increased'
		pause
	end if
	!write(*,*)' nzz=',nzz(ips)
	do izz=1,nzz
		read(1)((dv_3D(ips,ixx,iyy,izz),ixx=1,nxx),iyy=1,nyy)
	end do
	close(1)

	!write(*,*)' dv_3D(1,92,33,3)=',dv_3D(1,92,33,3)
	!pause

	!write(*,*)' nxx=',nxx,' nyy=',nyy,' nzz=',nzz
13	continue
end do

return
end