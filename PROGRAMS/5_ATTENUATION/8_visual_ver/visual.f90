USE DFPORT

character*4 dsaa/'DSAA'/
character*20 scale_line,scale_vpvs,scale,scale_vp,scale_vs
character*2 lv, ver
character*8 ar,md,line
character*1 itt,ps
allocatable dvan(:,:),vvv(:,:),vtmp(:,:)
real hlev(20)
real fia0(100),teta0(100),fib0(100),tetb0(100)
real fmark(200,50),tmark(200,50),smark(200,50)
integer kdot_rgb(3)
integer nmark(200)


common/pi/pi,per
common/center/fi0,tet0
common/keys/key_ft1_xy2

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0
rz=6371.


open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)ips
close(1)

write(ps,'(i1)')ips

write(*,*)' ar=',ar,' md=',md,' ps=',ps

i=system('mkdir ..\..\..\TMP_files\att')

!******************************************************************
key_ft1_xy2=1
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
read(1,*,end=441,err=441)key_true1      ! 1: true locations exist in separate line
read(1,*,end=441,err=441)key_flat1      ! 1: calculations in flat model, 2: spherical velocity
441 close(1)


!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=523)line
	if(line.eq.'ATTENUAT') goto 524
end do
523 continue
write(*,*)' cannot find ATTENUATION PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
524 read(1,*)
read(1,*) iter
read(1,*)val11,val12
read(1,*)val21,val22
att_aver=val11; add_att=val12
close(1)

if(ips.eq.2) att_aver=val21; add_att=val22
    


i=system('mkdir ..\..\..\TMP_files\vert')

key_preview=0
open(1,file='../../../preview_key.txt')
read(1,*,end=771)key_preview
771 close(1)

if(key_preview.ne.0) then
	i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\ATT')
557	open(1,file='../../../DATA/'//ar//'/config.txt')
	read(1,*) 
	read(1,*) npix_h_x0,npix_h_y0
	read(1,*) 
	read(1,*) 
	read(1,*) npix_x0,npix_y
	read(1,*)tick_x,tick_y
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*)scale_line
	read(1,*)a1,a2,dv_min0,dv_max0
	read(1,*)scale_vpvs
	read(1,*)vpvs_min,vpvs_max
	read(1,*,end=559,err=559)scale_vp
	read(1,*,end=559,err=559)vp_min,vp_max
	read(1,*,end=559,err=559)scale_vs
	read(1,*,end=559,err=559)vs_min,vs_max
	close(1)
    goto 558
559 write(*,*)'  '
    write(*,*)' Please update the file "config.txt" in the AREA folder !!!'
    write(*,*)' To the end of this file add lines with color scales for'
    write(*,*)' absolute P and S velocities. For example:'
    write(*,*)
    write(*,*)' rainbow_small.scl	scale for absolute Vp'
    write(*,*)' 2.6 3.2			diapason for absolute Vp'
    write(*,*)' rainbow_small.scl	scale for absolute Vs'
    write(*,*)' 1.7 2.3			diapason for absolute Vs'
    write(*,*)
    write(*,*)' When ready, press ENTER, and the program will make another attempt'
    close(1)
    pause
    goto 557
558 continue

	i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')
	i=system('copy ..\..\..\COMMON\scales_scl\'//scale_line//' '//scale_line)
	i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vpvs//' '//scale_vpvs)
	i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vp//' '//scale_vp)
	i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vs//' '//scale_vs)
end if

    
    
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
close(1)
!******************************************************************

ngr1=1
ngr2=nornt

!******************************************************************
if(key_ft1_xy2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
    do i=1,10000
        read(1,'(a8)',end=553)line
        if(line.eq.'AREA_CEN') goto 554
    end do
    553 continue
    write(*,*)' cannot find AREA CENTER in MAJOR_PARAM.DAT!!!'
    pause
    554 read(1,*)fi0,tet0
    close(1)
else
    fi0=0
    tet0=0
end if
!******************************************************************

open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=543)line
	if(line.eq.'ATTENUAT') goto 544
end do
543 continue
write(*,*)' cannot find ATTENUATION PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
544 read(1,*)
read(1,*)iter
close(1)

write(itt,'(i1)')iter

!******************************************************************


open(2,file='../../../DATA/'//ar//'/setver.dat')
read(2,*)nver
do ii=1,nver
	read(2,*) fia0(ii),teta0(ii),fib0(ii),tetb0(ii)
end do
read(2,*) dist_from_sec_event
read(2,*) dxsec
read(2,*) zmin,zmax,dzsec
read(2,*) dsmark
read(2,*) smaxx
read(2,*) ismth
551 close(2)

rsmth=ismth+0.5

call read_topo(ar)

do iver=1,nver
    write(ver,'(i2)')iver

    if(key_ft1_xy2.eq.1) then
        fia=fia0(iver)
        teta=teta0(iver)
        fib=fib0(iver)
        tetb=tetb0(iver)
        call SFDEC(fia,teta,0.,xa,ya,Z,fi0,tet0)
        call SFDEC(fib,tetb,0.,xb,yb,Z,fi0,tet0)
    else
        xa=fia0(iver)
        ya=teta0(iver)
        xb=fib0(iver)
        yb=tetb0(iver)
    end if
    write(*,*)' xa=',xa,' ya=',ya
    write(*,*)' xb=',xb,' yb=',yb
    dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
    !write(*,*)' dist=',dist

    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist
    nxsec=dist/dxsec+1
    dxsec=dist/(nxsec-1)
    nzsec=(zmax-zmin)/dzsec+1
    dzsec=(zmax-zmin)/(nzsec-1)
    write(*,*)' dist=',dist,' nxsec=',nxsec,' nzsec=',nzsec

    allocate (dvan(nxsec,nzsec),vvv(nxsec,nzsec),vtmp(nxsec,nzsec))
    vvv=0
    dvan=0

    open(11,file='../../../TMP_files/vert/mark_'//ver//'.dat')
    imark=0
    do sss=0.,dist,dsmark
        x=xa+cospov*sss
        y=ya+sinpov*sss
        if(key_ft1_xy2.eq.1) then
            call decsf(x,y,0.,fi0,tet0,FI,TET,h)
        else
            fi=x
            tet=y
        end if
        write(11,*)fi,tet,sss
        imark=imark+1
        fmark(iver,imark)=fi
        tmark(iver,imark)=tet
        smark(iver,imark)=sss
    end do
    imark=imark+1
    fmark(iver,imark)=fib
    tmark(iver,imark)=tetb
    smark(iver,imark)=dist
    close(11)
    nmark(iver)=imark

	
    ! Draw the position of the section on the surface (line)
    open(11,file='../../../TMP_files/vert/mark_'//ver//'.bln')
    write(11,*) imark
    do i=1,imark
	    write(11,*)fmark(iver,i),tmark(iver,i)	!,smark(i)
    end do
    close(11)

    !Draw topography on the section
    open(11,file='../../../TMP_files/vert/topo_'//ver//'.bln')
    write(11,*)nxsec
    do ix=1,nxsec
        sss=(ix-1)*dxsec
        xcur=xa+((xb-xa)/dist)*sss
        ycur=ya+((yb-ya)/dist)*sss
        if(key_ft1_xy2.eq.1) then
            call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
        else
            fff=xcur
            ttt=ycur
        end if
        depth=h_lim(fff,ttt)
        ztopo=flat_sph(key_flat1,xcur,ycur,depth)
        write(11,*)sss,-ztopo
    end do
    close(11)

    ! Read the coordinates of the stations
    open(2,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
    open(12,file='../../../TMP_files/vert/stat_'//ver//'.dat')
    i=0
    nst1=0
    3   i=i+1
        read(2,*,end=4)xst,yst,depth
        zst = flat_sph(key_flat1,xst,yst,depth)

        xx1=(xst-xa)*cospov+(yst-ya)*sinpov
        yy1=-(xst-xa)*sinpov+(yst-ya)*cospov

        if(abs(yy1).lt.dist_from_sec_event) then
            nst1=nst1+1
            write(12,*)xx1,-zst
        end if
        goto 3
    4   close(2)
    close(12)

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/srces'//itt//'.dat')
    open(11,file='../../../TMP_files/vert/ztr_'//ver//'.dat')
    nzt1=0
    nzt=0
872	read(1,*,end=871)fzt,tzt,depth
        if(key_ft1_xy2.eq.1) then
            call SFDEC(fzt,tzt,0.,xzt,yzt,Z,fi0,tet0)
        else
            xzt=fzt; yzt=tzt
        end if
        zzt = flat_sph(key_flat1,xzt,yzt,depth)
        nzt=nzt+1
        xx1=(xzt-xa)*cospov+(yzt-ya)*sinpov
        yy1=-(xzt-xa)*sinpov+(yzt-ya)*cospov
        if(abs(yy1).lt.dist_from_sec_event) then
            nzt1=nzt1+1
            write(11,*)xx1,-zzt,yy1
        end if
        goto 872
    871 close(1)
    close(11)
    write(*,*)' nst1=',nst1,' nzt1=',nzt1

    vvv=0
    dvan=0

!write(*,*)' dist=',dist,' nxsec=',nxsec

    do igr=ngr1,ngr2
        call prepare_model_att(ar,md,igr,ips)

        do ix=1,nxsec
            sss=(ix-1)*dxsec
            !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
            !sss=18
            !if(mod(ix,10).eq.0)write(*,*)' ix=',ix,' sss=',sss
            xcur=xa+((xb-xa)/dist)*sss
            ycur=ya+((yb-ya)/dist)*sss
   
            if(key_ft1_xy2.eq.1) then
                call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
            else
                fff=xcur; ttt=ycur
            end if
            !write(*,*)' xcur=',xcur,' ycur=',ycur
            !write(*,*)' fi=',fff,' tet=',ttt
            do iz=1,nzsec
                zcur=zmin+(iz-1)*dzsec
                !zcur=5
                depth = sph_flat(key_flat1,xcur,ycur,zcur)

                call att_1_grid(fff,ttt,zcur,smaxx, att,umn)
					
                depth = h_lim(fff,ttt)
                zlim_up = flat_sph(key_flat1,xcur,ycur,depth)



                if (zcur.lt.zlim_up) umn=0

                !write(*,*)' zcur=',zcur,' att=',att,' umn=',umn

                dvan(ix,iz) = dvan(ix,iz) + att*umn
                vvv(ix,iz)=vvv(ix,iz)+umn
            end do
        end do
    end do



    !***************************************************************
    !***************************************************************
    !***************************************************************

    dv_min=9999; dv_max=-9999
    do iz=nzsec,1,-1
        do ix=1,nxsec
            vanom=-999
            if(vvv(ix,iz).gt.0.0001) then
	            van=dvan(ix,iz)/vvv(ix,iz)
                vanom=van*100./att_aver
                if(vanom.lt.dv_min) dv_min=vanom
                if(vanom.gt.dv_max) dv_max=vanom
                !write(*,*)dv_min,dv_max
            end if
            dvan(ix,iz)=vanom
        end do
        !write(*,*)dv_min,dv_max
    end do
write(*,*)' dv_min=',dv_min,' dv_max=',dv_max,' att_aver=',att_aver

! smoothing:
        vtmp=dvan
        do iz=1,nzsec
	        do ix=1,nxsec
		        dvan(ix,iz)=-999
		        if(vvv(ix,iz).lt.0.01) cycle
		        vanm=0.
		        iv=0
		        do ixx=-ismth,ismth
			        if (ix+ixx.lt.1) cycle
			        if (ix+ixx.gt.nxsec) cycle
			        do izz=-ismth,ismth
				        if (iz+izz.lt.1) cycle
				        if (iz+izz.gt.nzsec) cycle
				        if(vvv(ix+ixx,iz+izz).lt.0.01) cycle
				        rr=ixx*ixx+izz*izz
				        r=sqrt(rr)
				        if(r.gt.rsmth) cycle
				        iv=iv+1
				        vanm=vanm+vtmp(ix+ixx,iz+izz)
			        end do
		        end do
		        dvan(ix,iz)=vanm/iv
	        end do
        end do


    open(14,file='../../../TMP_files/att/att_ver'//ps//ver//'.grd')
    write(14,'(a4)')dsaa
    write(14,*)nxsec,nzsec
    write(14,*)0,dist
    write(14,*)-zmax,-zmin
    write(14,*)dv_min,dv_max
    do iz=nzsec,1,-1
	    write(14,*)(dvan(ix,iz),ix=1,nxsec)
    end do
    close(14)

    deallocate (dvan,vvv,vtmp)

!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
        if(npix_x0.ne.0) then
	        npix_x=npix_x0
        else
	        npix_x = int(npix_y * dist/(zmax-zmin))
        end if
        write(*,*)' section:',ver,' dist=',dist,' npix_x=',npix_x
    
    
    
! PLOT THE ATTENUATION ANOMALIES
        open(14,file='config.txt')
        write(14,*)npix_x,npix_y
        write(14,*)'_______ Size of the picture in pixels (nx,ny)'
        write(14,*)0,dist
        write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
        write(14,*)-zmax,-zmin
        write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
        write(14,*)tick_x,tick_y
        write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
        write(14,51)ar,md,ps,ver
        51 format('..\..\..\PICS\',a8,'\',a8,'\ATT\ver_att',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        if(ips.eq.1) then
	        write(14,3431)	iver, iver
        3431 format(' P attenuation, section=',i1,'A - ',i1,'B')
        else
	        write(14,3432)	iver,iver
        3432 format(' S attenuation, section=',i1,'A - ',i1,'B')
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'

        59 format('********************************************')

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,52)ps,ver
        52 format('..\..\..\TMP_files\att\att_ver',a1,a2,'.grd')
        write(14,*)'_______ Location of the GRD file'
        write(14,53)scale_line
        53 format(a20)
        write(14,*)'_______ Scale for visualization'
        write(14,*)	dv_min0,dv_max0
        write(14,*)'_______ scale diapason'

        write(14,59)
        write(14,*)	2
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,2143)ver
        2143 format('..\..\..\TMP_files\vert\topo_',a2,'.bln')
        write(14,*)'_______ Location of the BLN file'
        write(14,*)	3
        write(14,*)'_______ Thickness of line in pixels'
        write(14,*)	0,0,0
        write(14,*)'_______ RGB color'

        close(14)

        i=system('layers.exe')


end do

stop
end