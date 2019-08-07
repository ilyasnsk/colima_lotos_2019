character*8 ar,md,line
character*1 it,gr,ps
real ornt(20)
real xray4(10000),yray4(10000),zray4(10000)
real plotray(:,:,:),dsy(:)

allocatable plotray,dsy

one=1.d0
pi=asin(one)*2.d0
per=pi/180.d0

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)ips		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(it,'(i1)')iter
write(gr,'(i1)')igr
write(ps,'(i1)')ips

write(*,*)' Computing ray density'
write(*,*)' ar=',ar,' md=',md,' ps=',ps,' gr=',gr

!******************************************************************
key_ft1_xy2=1
key_true1=0
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
    read(1,'(a8)',end=513)line
    if(line.eq.'GENERAL ') goto 514
end do
513 continue
write(*,*)' cannot find GENERAL INFORMATION in MAJOR_PARAM.DAT!!!'
pause
514 continue
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*,end=441,err=441)key_ft1_xy2
read(1,*,end=441,err=441)key_true1
441 close(1)
!******************************************************************

open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=573)line
	if(line.eq.'ORIENTAT') goto 574
end do
573 continue
write(*,*)' cannot find ORIENTATIONS in MAJOR_PARAM.DAT!!!'
pause
574 read(1,*)nornt
read(1,*)(ornt(i),i=1,nornt)
close(1)

!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=533)line
	if(line.eq.'GRID_PAR') goto 534
end do
533 continue
write(*,*)' cannot find GRID_PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
534 continue
read(1,*)xlim1,xlim2,dxpl
read(1,*)ylim1,ylim2,dypl
read(1,*)zlim1,zlim2,dzpl
close(1)
!******************************************************************


orient=ornt(igr)
write(*,*)' orient=',orient
sinbase=sin(orient*per)
cosbase=cos(orient*per)

nxpl=int_best((xlim2-xlim1)/dxpl)
nypl=int_best((ylim2-ylim1)/dypl)
nzpl=int_best((zlim2-zlim1)/dzpl)

write(*,*)' nx=',nxpl,' ny=',nypl,' nz=',nzpl

allocate(plotray(nxpl,nypl,nzpl),dsy(nypl))

plotray=0.
dsy=0.
nray=0


open(1,file='../../../DATA/'//ar//'/'//md//'/data/atten'//ps//'.dat')
open(2,file='../../../DATA/'//ar//'/'//md//'/data/ray_paths_att'//ps//'.dat',form='binary')

open(31,file='../../../DATA/'//ar//'/'//md//'/data/dist_atten'//ps//'.dat')


nray=0
aver=0

728 continue
    read(1,*,end=729)xzt,yzt,zzt,xst,yst,zst,att,tout
    nray=nray+1

    read(2)npray
    !write(*,*)npray
    if(npray.eq.0)goto 728
    do i=1,npray
	    read(2)xray4(i),yray4(i),zray4(i)
    end do
!if(nzt.lt.5814) cycle
!write(*,*)ikr,ips,npray,nkrat
!if(ikr.lt.71) cycle
        
    dist=0
    do ip=2,npray
        xl1=xray4(ip-1)
        xl2=xray4(ip)
        yl1=yray4(ip-1)
        yl2=yray4(ip)
        zl1=zray4(ip-1)
        zl2=zray4(ip)
        dss=sqrt((xl1-xl2)*(xl1-xl2)+(yl1-yl2)*(yl1-yl2)+(zl1-zl2)*(zl1-zl2))
        dist=dist+dss

        xmid=(xl1+xl2)/2.
        ymid=(yl1+yl2)/2.
        zmid=(zl1+zl2)/2.
        !write(*,*)ip,xmid,ymid,zmid

        xxx=xmid*cosbase+ymid*sinbase
        yyy=-xmid*sinbase+ymid*cosbase
        xmid=xxx
        ymid=yyy

        if((xmid-xlim1)*(xmid-xlim2).gt.0) cycle
        if((ymid-(ylim1+0.5*dypl))*(ymid-(ylim2-0.5*dypl)).gt.0) cycle
        if((zmid-zlim1)*(zmid-zlim2).gt.0) cycle

!xmid=0
!ymid=0
!zmid=60

        ilevy=0
        do iy=1,nypl
	        y1=ylim1+(iy-0.5)*dypl
	        y2=ylim1+(iy+0.5)*dypl
	        !write(*,*)' y1=',y1,' y2=',y2
	        if((ymid-y1)*(ymid-y2).le.0.)exit
        end do

        if(iy.eq.0.or.iy.gt.nypl) cycle

        do ix=1,nxpl
	        x1=xlim1+dxpl*(ix-1)
	        x2=xlim1+dxpl*ix
	        if((xmid-x1)*(xmid-x2).le.0.)exit
        end do

        if(ix.eq.0.or.ix.gt.nxpl) cycle

        do iz=1,nzpl
	        z1=zlim1+dzpl*(iz-1)
	        z2=zlim1+dzpl*iz
	        if((zmid-z1)*(zmid-z2).le.0.)exit
        end do

        if(iz.eq.0.or.iz.gt.nzpl) cycle



        plotray(ix,iy,iz)=plotray(ix,iy,iz)+dss
        dsy(iy)=dsy(iy)+dss

    end do
    write(31,*)dist,att,10./att
    aver=aver+att/dist

    if(mod(nray,100).eq.0) write(*,*)' nray=',nray,' dist=',dist,' att=',att

    goto 728
729 continue
close(1)
close(2)
close(31)
write(*,*)' total number of rays:',nray




open(15,file='../../../DATA/'//ar//'/'//md//'/data/num_att'//ps//'.dat')
write(15,*) nray
close(15)

open(12,file='../../../TMP_files/tmp/plot_att'//gr//ps//'.dat')
do iy=1,nypl
	write(12,*)((plotray(ix,iy,iz),ix=1,nxpl),iz=1,nzpl)
end do
close(12)
stop
end

