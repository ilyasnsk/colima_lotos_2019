subroutine read_topo(ar)
character*8 md,ar
common/topo/nxx,xmin,xmax,nyy,ymin,ymax,htopo(1500,1500),topomin,topomax

nfi=0
topomin=0
topomax=0

open(1,file='../../DATA/'//ar//'/inidata/topo_xyz.grd',status='old',err=1)

read(1,*)
read(1,*)nxx,nyy
read(1,*)xmin,xmax
read(1,*)ymin,ymax
read(1,*)hmin,hmax
do iyy=1,nyy
	read(1,*)(htopo(ixx,iyy),ixx=1,nxx)
end do
close(1)

topomin=999999
topomax=-999999

do ixx=1,nxx
    do iyy=1,nyy
        if(htopo(ixx,iyy).gt.topomax) topomax=htopo(ixx,iyy)
        if(htopo(ixx,iyy).lt.topomin) topomin=htopo(ixx,iyy)
    end do
end do
!write(*,*)' nar=',nar

1 continue

!write(*,*)' nxx=',nxx,' nyy=',nyy
!write(*,*)' xmin=',xmin,' xmax=',xmax
!write(*,*)' ymin=',ymin,' ymax=',ymax
!write(*,*)' topomin=',topomin,' topomax=',topomax
return
end