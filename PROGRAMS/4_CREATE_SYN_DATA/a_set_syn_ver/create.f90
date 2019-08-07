! Visualization of synthetic model in VERTICAL sections

USE DFPORT

character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*20 scale_dv,scale_vpvs,scale,scale_vp,scale_vs
character*1 ps
character*2 ver,lv
real fia0(100),teta0(100),fib0(100),tetb0(100)
real hlev(100)
real fmark(200,20),tmark(200,20),smark(200,20)
integer nmark(100)

integer kdot_rgb(3)
common/center/fi0,tet0
common/pi/pi,per
common/keys/key_ft1_xy2

allocatable v_ini(:,:),v_abs_ini(:,:,:)


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

i=system('mkdir ..\..\..\TMP_files\vert')


open(1,file='../../../model.dat')
read(1,'(a8)')ar	! synthetic model
read(1,'(a8)')md	! synthetic model
close(1)

write(*,*)' ar=',ar,' md=',md

key_preview=0
open(1,file='../../../preview_key.txt')
read(1,*,end=771)key_preview
771 close(1)

if(key_preview.ne.0) then
    i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\SYN')

557 open(1,file='../../../DATA/'//ar//'/config.txt')
    read(1,*) 
    read(1,*) npix_h_x0,npix_h_y0
    read(1,*) tick_h_x,tick_h_y
    read(1,*) 
    read(1,*) npix_x0,npix_y
    read(1,*)tick_x,tick_y
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)scale_dv
    read(1,*)dv_min,dv_max
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
    i=system('copy ..\..\..\COMMON\scales_scl\'//scale_dv//' '//scale_dv)
    i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vpvs//' '//scale_vpvs)
    i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vp//' '//scale_vp)
    i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vs//' '//scale_vs)

end if


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
read(1,*,end=741,err=741)
read(1,*,end=741,err=741)
read(1,*,end=741,err=741)key_ft1_xy2
741 close(1)
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

call read_topo(ar)
call read_vref(ar,md)
call read_anom(ar,md)


ksec_ver=1
kps_ver=1
open(2,file='../../../DATA/'//ar//'/setver.dat')
read(2,*)nver
do ii=1,nver
    read(2,*) fia0(ii),teta0(ii),fib0(ii),tetb0(ii)
end do
read(2,*) 
read(2,*) dxsec
read(2,*) zmin,zmax,dzsec
read(2,*) dsmark
read(2,*) dismax
read(2,*) ismth
read(2,*) 
read(2,*,end=551) ksec_ver,kps_ver
read(2,*,end=551) dfmark,dtmark
close(2)
!dxsec=2
!dzsec=2

551 close(2)


open(16,file='../../../TMP_files/info_map.txt')
write(16,*)nver
write(16,*)fmap1+fi0,fmap2+fi0
write(16,*)tmap1+tet0,tmap2+tet0


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
    !write(*,*)' xa=',xa,' ya=',ya
    !write(*,*)' xb=',xb,' yb=',yb
    dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
    npix_x = int(npix_y * dist/(zmax-zmin))
    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist
    nxsec=dist/dxsec+1
    dxsec=dist/(nxsec-1)
    nzsec=(zmax-zmin)/dzsec+1
    dzsec=(zmax-zmin)/(nzsec-1)


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
    if(key_ft1_xy2.eq.1) then
        fmark(iver,imark)=fib
        tmark(iver,imark)=tetb
    else
        fmark(iver,imark)=xb
        tmark(iver,imark)=yb
    end if
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

    write(16,*)'mark_',ver,'.bln'

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
        topo=h_lim(fff,ttt)
        write(11,*)sss,-topo
    end do
    close(11)


    allocate (v_ini(nxsec,nzsec),v_abs_ini(2,nxsec,nzsec))

    write(*,*)' section:',ver,' dist=',dist,' nxsec=',nxsec,' nzsec=',nzsec
    v_ini=0
    v_abs_ini=0
	
    do ips=1,2
        write(ps,'(i1)')ips

        do ix=1,nxsec
            sss=(ix-1)*dxsec
            !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
            !sss=40
            !if(mod(ix,10).eq.0)write(*,*)' ix=',ix,' sss=',sss
            xcur=xa+((xb-xa)/dist)*sss
            ycur=ya+((yb-ya)/dist)*sss
            if(key_ft1_xy2.eq.1) then
                call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
            else
                fff=xcur
                ttt=ycur
            end if

            relief=h_lim(fff,ttt)



            !write(*,*)' xcur=',xcur,' ycur=',ycur
            !write(*,*)' fi=',fi,' tet=',tet,' relief=',relief
            do iz=1,nzsec
                zcur=zmin+iz*dzsec
                v0=vrefmod(zcur,ips)
                dv=anomaly(xcur,ycur,zcur,ips)

                if (zcur.lt.relief) then
                    v_ini(ix,iz)=-999
                    v_abs_ini(ips,ix,iz)=-999
                else
                    v_ini(ix,iz)=dv
                    v_abs_ini(ips,ix,iz)=v0*(1+dv/100)
                end if
            end do

        end do

       open(14,file='../../../TMP_files/vert/syn_dv'//ps//ver//'.xyz')
       open(15,file='../../../TMP_files/vert/syn_abs'//ps//ver//'.xyz')
       do ix=1,nxsec
            sss=(ix-1)*dxsec
            do iz=1,nzsec
                zcur=zmin+iz*dzsec
                write(14,*)sss,-zcur,v_ini(ix,iz)
                write(15,*)sss,-zcur,v_abs_ini(ips,ix,iz)
           end do
        end do
        close(14)
        close(15)




        open(14,file='../../../TMP_files/vert/syn_dv'//ps//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
            write(14,*)(v_ini(ix,iz),ix=1,nxsec)
        end do
        close(14)

        open(14,file='../../../TMP_files/vert/syn_abs'//ps//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
            write(14,*)(v_abs_ini(ips,ix,iz),ix=1,nxsec)
        end do
        close(14)


        if(key_preview.eq.0) goto 441


        !*********************************************************
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
        !write(*,*)' section:',ver,' dist=',dist,' npix_x=',npix_x

!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
! PLOT THE VELOCITY ANOMALIES

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
        51 format('..\..\..\PICS\',a8,'\',a8,'\SYN\vert_syn_dv',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        if (ips.eq.1) then
            write(14,3431)	iver, iver
            3431 format(' Synthetic P anomalies, section=',i1,'A - ',i1,'B')
        else
            write(14,3432)	iver, iver
            3432 format(' Synthetic S anomalies, section=',i1,'A - ',i1,'B')
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'

        59 format('********************************************')

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,52)ps,ver
        52 format('..\..\..\TMP_files\vert\syn_dv',a1,a2,'.grd')
        write(14,*)'_______ Location of the GRD file'
        write(14,53)scale_dv
        53 format(a20)
        write(14,*)'_______ Scale for visualization'
        write(14,*)	dv_min,dv_max
        write(14,*)'_______ scale diapason'

        close(14)


        i=system('layers.exe')


!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
! PLOT THE ABSOLUTE VELOCITIES

        amp_min=vp_min
        amp_max=vp_max
        scale=scale_vp
        if(ips.eq.2) then
            amp_min=vs_min
            amp_max=vs_max
            scale=scale_vs
        end if


        open(14,file='config.txt')
        write(14,*)npix_x,npix_y
        write(14,*)'_______ Size of the picture in pixels (nx,ny)'
        write(14,*)0,dist
        write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
        write(14,*)-zmax,-zmin
        write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
        write(14,*)tick_x,tick_y
        write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
        write(14,517)ar,md,ps,ver
        517 format('..\..\..\PICS\',a8,'\',a8,'\SYN\vert_syn_abs',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        if (ips.eq.1) then
            write(14,347)	iver, iver
            347 format(' Synthetic P absolute, section=',i1,'A - ',i1,'B')
        else
            write(14,732)	iver, iver
            732 format(' Synthetic S absolute, section=',i1,'A - ',i1,'B')
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,752)ps,ver
        752 format('..\..\..\TMP_files\vert\syn_abs',a1,a2,'.grd')
        write(14,*)'_______ Location of the GRD file'
        write(14,753)scale
        753 format(a20)
        write(14,*)'_______ Scale for visualization'
        write(14,*)	amp_min,amp_max
        write(14,*)'_______ scale diapason'

        close(14)

        i=system('layers.exe')

441	continue

    end do      ! ips=1,2

    v_ini=0

    !write(*,*)' only  S:',v_abs_ini(1,135,183),v_abs_ini(2,135,183),v_abs_ini(1,135,183)/v_abs_ini(2,135,183)
    !write(*,*)' P and S:',v_abs_ini(1,151,133),v_abs_ini(2,151,133),v_abs_ini(1,151,133)/v_abs_ini(2,151,133)
    !pause


    do ix=1,nxsec
        sss=(ix-1)*dxsec

        do iz=1,nzsec
            zcur=zmin+iz*dzsec
            !write(*,*)sss,zcur,v_abs_ini(1,ix,iz),v_abs_ini(2,ix,iz),v_abs_ini(1,ix,iz)/v_abs_ini(2,ix,iz)
            v_ini(ix,iz)=v_abs_ini(1,ix,iz)/v_abs_ini(2,ix,iz)
        end do
        !pause
    end do

    open(14,file='../../../TMP_files/vert/syn_vpvs'//ver//'.xyz')
    do ix=1,nxsec
        sss=(ix-1)*dxsec
        do iz=1,nzsec
            zcur=zmin+iz*dzsec
            write(14,*)sss,-zcur,v_ini(ix,iz)
        end do
    end do
    close(14)

    open(14,file='../../../TMP_files/vert/syn_vpvs'//ver//'.grd')
    write(14,'(a4)')dsaa
    write(14,*)nxsec,nzsec
    write(14,*)0,dist
    write(14,*)-zmax,-zmin
    write(14,*)-9999,999
    do iz=nzsec,1,-1
        write(14,*)(v_ini(ix,iz),ix=1,nxsec)
    end do
    close(14)


    if(key_preview.eq.0) goto 442


    !*********************************************************
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
    !write(*,*)' section:',ver,' dist=',dist,' npix_x=',npix_x

    open(14,file='config.txt')
    write(14,*)npix_x,npix_y
    write(14,*)'_______ Size of the picture in pixels (nx,ny)'
    write(14,*)0,dist
    write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
    write(14,*)-zmax,-zmin
    write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
    write(14,*)tick_x,tick_y
    write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
    write(14,61)ar,md,ver
    61 format('..\..\..\PICS\',a8,'\',a8,'\SYN\vert_syn_vpvs',a2,'.png')
    write(14,*)'_______ Path of the output picture'
    write(14,3433)	iver, iver
    3433 format(' Synthetic Vp/Vs ratio, section=',i1,'A - ',i1,'B')
    write(14,*)'_______ Title of the plot on the upper axe'
    write(14,*)	1
    write(14,*)'_______ Number of layers'

    69 format('********************************************')

    write(14,69)
    write(14,*)	1
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,62)ver
    62 format('..\..\..\TMP_files\vert\syn_vpvs',a2,'.grd')
    write(14,*)'_______ Location of the GRD file'
    write(14,63)scale_vpvs
    63 format(a20)
    write(14,*)'_______ Scale for visualization'
    write(14,*)	vpvs_min,vpvs_max
    write(14,*)'_______ scale diapason'


    close(14)


    i=system('layers.exe')


    442		continue

    deallocate(v_ini,v_abs_ini)

end do      ! iver=1,nver
close(16)


! LOCATIONS OF THE PROFILES IN MAP VIEW:

if(key_preview.eq.0) stop

open(2,file='../../../DATA/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
close(2)

cost=cos(tet0*per)
if(npix_h_x0.ne.0) then
    npix_x=npix_h_x0
else
    npix_x=int(npix_h_y0*((fmap2-fmap1)/(tmap2-tmap1))*cost)
    npix_y=npix_h_y0
end if

if(npix_h_y0.ne.0) then
    npix_y=npix_h_y0
else
    npix_y=int(npix_h_x0*((tmap2-tmap1)/(fmap2-fmap1))/cost)
    npix_x=npix_h_x0
end if

write(*,*)' npix_x=',npix_x,' npix_y=',npix_y



izmap=ksec_ver
izzz=hlev(izmap)
iiips=kps_ver

write(ps,'(i1)')iiips
write(lv,'(i2)')izmap

open(14,file='config.txt')
write(14,*)npix_x,npix_y
write(14,*)'_______ Size of the picture in pixels (nx,ny)'
write(14,*)fmap1+fi0,fmap2+fi0
write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
write(14,*)tmap1+tet0,tmap2+tet0
write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
write(14,*)tick_h_x,tick_h_y
write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
write(14,151)ar,md,ps,lv
151 format('..\..\..\PICS\',a8,'\',a8,'\SYN\hor_profiles',a1,a2,'.png')
write(14,*)'_______ Path of the output picture'
if(iiips.eq.1) then
	write(14,5331)	izzz
5331 format(' P anomalies, depth=',i4)
else
	write(14,5332)	izzz
5332 format(' S anomalies, depth=',i4)
end if
write(14,*)'_______ Title of the plot on the upper axe'
write(14,*)	1
write(14,*)'_______ Number of layers'

write(14,59)
write(14,*)	1
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,552)ps,lv
552 format('..\..\..\TMP_files\hor\syn',a1,'_',a2,'.grd')
write(14,*)'_______ Location of the GRD file'
write(14,555)scale_dv
555 format(a20)
write(14,*)'_______ Scale for visualization'
write(14,*)	dv_min,dv_max
write(14,*)'_______ scale diapason'

open(1,file='../../../DATA/'//ar//'/map/polit_bound.bln',status='old',err=491)
close(1)
write(14,59)
write(14,*)	2
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,54)ar
54 format('..\..\..\DATA\',a8,'\map\polit_bound.bln')
write(14,*)'_______ Location of the BLN file'
write(14,*)	3
write(14,*)'_______ Thickness of line in pixels'
write(14,*)	100,0,0
write(14,*)'_______ RGB color'
491 continue

open(1,file='../../../DATA/'//ar//'/map/coastal_line.bln',status='old',err=492)
close(1)
write(14,59)
write(14,*)	2
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,55)ar
55 format('..\..\..\DATA\',a8,'\map\coastal_line.bln')
write(14,*)'_______ Location of the BLN file'
write(14,*)	3
write(14,*)'_______ Thickness of line in pixels'
write(14,*)	0,0,0
write(14,*)'_______ RGB color'
492 continue

do iver=1,nver
	write(ver,'(i2)')iver

	write(14,59)
	write(14,*)	2
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,154)ver
154 format('..\..\..\TMP_files\vert\mark_',a2,'.bln')
	write(14,*)'_______ Location of the BLN file'
	write(14,*)	2
	write(14,*)'_______ Thickness of line in pixels'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	write(14,59)
	write(14,*)	3
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,156)ver
156 format('..\..\..\TMP_files\vert\mark_',a2,'.dat')
	write(14,*)'_______ Location of the DAT file'
	write(14,*)	1
	write(14,*)'_______ Symbol (1: circle, 2: square)'
	write(14,*)	6
	write(14,*)'_______ Size of dots in pixels'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	fia=fmark(iver,1)+dfmark
	teta=tmark(iver,1)+dtmark

	write(14,59)
	write(14,*)	5
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,756)iver
756 format(i2,' - A')
	write(14,*)'_______ Text to output'
	write(14,*)	fia,teta
	write(14,*)'_______ Coordinates'
	write(14,*)	1
	write(14,*)'_______ Type (1: Arial, 2: Times New Romans)'
	write(14,*)	12
	write(14,*)'_______ Size'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	fib=fmark(iver,nmark(iver))+dfmark
	tetb=tmark(iver,nmark(iver))+dtmark

	write(14,59)
	write(14,*)	5
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,757)iver
757 format(i2,' - B')
	write(14,*)'_______ Text to output'
	write(14,*)	fib,tetb
	write(14,*)'_______ Coordinates'
	write(14,*)	1
	write(14,*)'_______ Type (1: Arial, 2: Times New Romans)'
	write(14,*)	12
	write(14,*)'_______ Size'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	do imark=2,nmark(iver)-2

		fi=fmark(iver,imark)+dfmark
		tet=tmark(iver,imark)+dtmark
		isss=smark(iver,imark)

		write(14,59)
		write(14,*)	5
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,'(i4)')isss
		write(14,*)'_______ Text to output'
		write(14,*)	fi,tet
		write(14,*)'_______ Coordinates'
		write(14,*)	1
		write(14,*)'_______ Type (1: Arial, 2: Times New Romans)'
		write(14,*)	10
		write(14,*)'_______ Size'
		write(14,*)	0,0,0
		write(14,*)'_______ RGB color'

	end do

end do

close(14)


i=system('layers.exe')


stop
end