USE DFPORT

character*4 dsaa/'DSAA'/
character*20 scale_line,scale_vpvs,scale,scale_vp,scale_vs

character*8 ar,md,line
character*2 lv, ver
character*1 ps ,rm,it
allocatable dvan(:,:),aprio_dv(:,:),vvv(:,:),vtmp(:,:),vabs(:,:),vabsp(:,:)
real hlev(20)
real fia0(100),teta0(100),fib0(100),tetb0(100)
real fmark(200,20),tmark(200,20),smark(200,20)
integer nmark(100)


common/pi/pi,per
common/center/fi0,tet0
common/keys/key_ft1_xy2

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0
rz=6371.

w_limit=0.2

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)iter		
close(1)

!md='MODEL_11'; iter=5
write(it,'(i1)')iter

write(*,*)
write(*,*)' ***********************************************'
write(*,*)' VISUALISATION in vertical section: '
write(*,*)' ar=',ar,' md=',md,' iter=',iter

!******************************************************************
key_ft1_xy2=1
key_true1=0
key_flat1=1
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
read(1,*,end=441,err=441)key_ft1_xy2    ! 1: geogr, 2: Cartesian, 
read(1,*,end=441,err=441)key_true1      ! 1: true locations exist in separate line
read(1,*,end=441,err=441)key_flat1      ! 1: calculations in flat model, 2: spherical velocity
441 close(1)
!******************************************************************


i=system('mkdir ..\..\..\TMP_files\vert')

key_preview=0
open(1,file='../../../preview_key.txt')
read(1,*,end=771)key_preview
771 close(1)

if(key_preview.ne.0) then
	i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\IT'//it)
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


add_perc=0
kod_av_bias=0
kod_apriori=0
ind_srce=1


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

ksec_ver=1
kps_ver=1

open(2,file='../../../DATA/'//ar//'/setver.dat')
read(2,*)nver
do ii=1,nver
	read(2,*) fia0(ii),teta0(ii),fib0(ii),tetb0(ii)
end do
read(2,*) dist_from_sec_event
read(2,*) dxsec
read(2,*) zmin,zmax,dzsec
read(2,*) dsmark
read(2,*) dismax
read(2,*) ismth
read(2,*) 
read(2,*,end=551) ksec_ver,kps_ver
read(2,*,end=551) dfmark,dtmark
551 close(2)


write(it,'(i1)')iter


rsmth=ismth+0.5


call read_vref(ar,md)


if(kod_apriori.eq.1) then
	call read_ini_model(ar,md)
end if

call read_topo(ar)

open(2,file='../../../DATA/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
close(2)

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

    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist
    nxsec=dist/dxsec+1
    dxsec=dist/(nxsec-1)
    nzsec=(zmax-zmin)/dzsec+1
    dzsec=(zmax-zmin)/(nzsec-1)
    !write(*,*)' dist=',dist,' nxsec=',nxsec

    allocate (dvan(nxsec,nzsec),aprio_dv(nxsec,nzsec),vvv(nxsec,nzsec))
    allocate (vtmp(nxsec,nzsec),vabs(nxsec,nzsec),vabsp(nxsec,nzsec))
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
        depth=h_lim(fff,ttt)
        ztopo=flat_sph(key_flat1,xcur,ycur,depth)
        write(11,*)sss,-ztopo
    end do
    close(11)

    if(ind_srce.ne.0) then
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



        !open(1,file='../../data/'//ar//'/inidata/events.dat')

        open(1,file='../../../DATA/'//ar//'/'//md//'/data/srces'//it//'.dat')
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
    end if

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
    do ips=1,2
        write(ps,'(i1)')ips
        vvv=0
        dvan=0

    !write(*,*)' dist=',dist,' nxsec=',nxsec

        do igr=ngr1,ngr2
            call prepare_model_v(ar,md,ips,iter,igr)

            do ix=1,nxsec
                sss=(ix-1)*dxsec
                !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
                !sss=500
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
                    depth = sph_flat(key_flat1,xcur,ycur,zcur)
                    vref=vrefmod(depth,ips)
                    if(igr.eq.ngr1) then
	                    aprio_dv(ix,iz)=0
	                    if(kod_apriori.eq.1) then
		                    aprio_dv(ix,iz)=vert_anom(xcur,ycur,zcur,ips)
	                    end if
                    end if
                    !zcur=15
                    !write(*,*)' avmoho=',avmoho,' z_moho=',z_moho

                    call dv_1_grid_v(fff,ttt,zcur,dismax,   dv,umn)
					
                    depth = h_lim(fff,ttt)
                    zlim_up = flat_sph(key_flat1,xcur,ycur,depth)



                    if (zcur.lt.zlim_up) umn=0

                    !write(*,*)' zcur=',zcur,' dv=',dv,' umn=',umn
                    !pause

                    dvan(ix,iz) = dvan(ix,iz) + dv*umn
                    vvv(ix,iz)=vvv(ix,iz)+umn
                end do
            end do
        end do

!        write(*,*)' nxsec=',nxsec,' nzsec=',nzsec



        !***************************************************************
        !***************************************************************
        !***************************************************************

        do iz=nzsec,1,-1
        !write(*,*)' iz=',iz
            do ix=1,nxsec
                vanom=-999
                if(vvv(ix,iz).gt.0.0001) then
	                vanom=dvan(ix,iz)/vvv(ix,iz)
                end if
                dvan(ix,iz)=vanom
            end do
        end do

!        write(*,*)' 2: nxsec=',nxsec,' nzsec=',nzsec

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

!        write(*,*)' 3: nxsec=',nxsec,' nzsec=',nzsec


        aver=0
        naver=0
        do iz=1,nzsec
            zcur=zmin+iz*dzsec
            do ix=1,nxsec
                sss=(ix-1)*dxsec
                !write(*,*)' ix=',ix,' sss=',sss,' dist=',dist
                !sss=500
                !if(mod(ix,10).eq.0)write(*,*)' ix=',ix,' sss=',sss
                xcur=xa+((xb-xa)/dist)*sss
                ycur=ya+((yb-ya)/dist)*sss
                depth = sph_flat(key_flat1,xcur,ycur,zcur)
                v0=vrefmod(depth,ips)
                vanom=-999
                vab=-999
                if(vvv(ix,iz).gt.w_limit) then
                    vanom=100*dvan(ix,iz)/v0
                    dvan(ix,iz)=vanom
                    vab=v0*(1+0.01*(vanom+aprio_dv(ix,iz)))
                    !if(ips.eq.2) then
                    !write(*,*)' vanom=',vanom,' vab=',vab
                    !pause
                    !end if
                    aver=aver+vanom
                    naver=naver+1
                end if
                vabs(ix,iz)=vab
            end do
        end do
!        write(*,*)' 4: nxsec=',nxsec,' nzsec=',nzsec
                   !if(ips.eq.2) then
                        !write(*,*)((vabs(ix,iz),ix=1,nxsec),iz=1,nzsec)
                        !stop
                    !end if

        !pause
        aver=aver/naver
        !dvan=dvan-aver

        if(ips.eq.1) then
	        vabsp=vabs
        else
	        do iz=1,nzsec
		        !zcur=zmin+iz*dzsec
		        !write(*,*)' zcur=',zcur
		        do ix=1,nxsec
			        !sss=(ix-1)*dxsec
			        !write(*,*)' sss=',sss
			        !write(*,*)' vp=',vabsp(ix,iz),' vs=',vabs(ix,iz),' vp/vs=',vabsp(ix,iz)/vabs(ix,iz)
			        if(abs(vabs(ix,iz)).gt.900.or.abs(vabsp(ix,iz)).gt.900) then
				        vtmp(ix,iz)=-999
			        else
				        vtmp(ix,iz)=vabsp(ix,iz)/vabs(ix,iz)
				        !write(*,*)' vp=',vabsp(ix,iz),' vs=',vabs(ix,iz),' vp/vs=',vabsp(ix,iz)/vabs(ix,iz)
			        end if

			        !vabs(ix,iz)=vab
		        end do
	        end do


	        open(14,file='../../../TMP_files/vert/vpvs_'//it//ver//'.xyz')
	        do iz=nzsec,1,-1
		        zcur=zmin+iz*dzsec
		        do ix=1,nxsec
			      sss=(ix-1)*dxsec
                            write(14,*)sss,-zcur,vtmp(ix,iz)
		        end do
	        end do
               close(14)

    !write(*,*)' ****************** dist=',dist,' nxsec=',nxsec


	        open(14,file='../../../TMP_files/vert/vpvs_'//it//ver//'.grd')
	        write(14,'(a4)')dsaa
	        write(14,*)nxsec,nzsec
	        write(14,*)0,dist
	        write(14,*)-zmax,-zmin
	        write(14,*)-9999,999
	        do iz=nzsec,1,-1
		        write(14,*)(vtmp(ix,iz),ix=1,nxsec)
	        end do
	        close(14)

			
        end if

        open(14,file='../../../TMP_files/vert/anom_'//it//ver//'.xyz')
        open(15,file='../../../TMP_files/vert/abs_'//it//ver//'.xyz')
        do iz=nzsec,1,-1
	        zcur=zmin+iz*dzsec
	        do ix=1,nxsec
		        sss=(ix-1)*dxsec
                        write(14,*)sss,-zcur,dvan(ix,iz)+add_perc+aprio_dv(ix,iz)
                        write(15,*)sss,-zcur,vabs(ix,iz)
	        end do
        end do
        close(14)
        close(15)


        open(14,file='../../../TMP_files/vert/ver_'//ps//it//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
	        write(14,*)(dvan(ix,iz)+add_perc+aprio_dv(ix,iz),ix=1,nxsec)
        end do
        close(14)

        open(14,file='../../../TMP_files/vert/abs_'//ps//it//ver//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nxsec,nzsec
        write(14,*)0,dist
        write(14,*)-zmax,-zmin
        write(14,*)-9999,999
        do iz=nzsec,1,-1
	        write(14,*)(vabs(ix,iz),ix=1,nxsec)
        end do
        close(14)

        if(key_preview.eq.0) goto 541


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
        write(*,*)' section:',ver,' dist=',dist,' npix_x=',npix_x

!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
!*********************************************************
! PLOT THE VELOCITY ANOMALIES
        !open(14,file='../../CREATE_PICS/LAYERS/config.txt')
        open(14,file='config.txt')
        write(14,*)npix_x,npix_y
        write(14,*)'_______ Size of the picture in pixels (nx,ny)'
        write(14,*)0,dist
        write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
        write(14,*)-zmax,-zmin
        write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
        write(14,*)tick_x,tick_y
        write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
        write(14,51)ar,md,it,ps,ver
        51 format('..\..\..\PICS\',a8,'\',a8,'\IT',a1,'\ver_dv',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        if(ips.eq.1) then
	        write(14,3431)	iver, iver
        3431 format(' P anomalies, section=',i1,'A - ',i1,'B')
        else
	        write(14,3432)	iver,iver
        3432 format(' S anomalies, section=',i1,'A - ',i1,'B')
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'

        59 format('********************************************')

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,52)ps,it,ver
        52 format('..\..\..\TMP_files\vert\ver_',a1,a1,a2,'.grd')
        write(14,*)'_______ Location of the GRD file'
        write(14,53)scale_line
        53 format(a20)
        write(14,*)'_______ Scale for visualization'
        write(14,*)	dv_min,dv_max
        write(14,*)'_______ scale diapason'

        write(14,59)
        write(14,*)	3
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,56)ver
        56 format('..\..\..\TMP_files\vert\ztr_',a2,'.dat')
        write(14,*)'_______ Location of the DAT file'
        write(14,*)	1
        write(14,*)'_______ Symbol (1: circle, 2: square)'
        write(14,*)	5
        write(14,*)'_______ Size of dots in pixels'
        write(14,*)	250,0,0
        write(14,*)'_______ RGB color'

        write(14,59)
        write(14,*)	3
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,76)ver
        76 format('..\..\..\TMP_files\vert\stat_',a2,'.dat')
        write(14,*)'_______ Location of the DAT file'
        write(14,*)	2
        write(14,*)'_______ Symbol (1: circle, 2: square)'
        write(14,*)	7
        write(14,*)'_______ Size of dots in pixels'
        write(14,*)	0,0,250
        write(14,*)'_______ RGB color'


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
        write(14,510)ar,md,it,ps,ver
        510 format('..\..\..\PICS\',a8,'\',a8,'\IT',a1,'\abs_ver',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        if(ips.eq.1) then
	        write(14,340)	iver, iver
        340 format(' P absolute, section=',i1,'A - ',i1,'B')
        else
	        write(14,302)	iver,iver
        302 format(' S absolute, section=',i1,'A - ',i1,'B')
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'



        amp_min=vp_min
        amp_max=vp_max
        scale=scale_vp
        if(ips.eq.2) then
            amp_min=vs_min
            amp_max=vs_max
            scale=scale_vs
        end if

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,520)ps,it,ver
        520 format('..\..\..\TMP_files\vert\abs_',a1,a1,a2,'.grd')
        write(14,*)'_______ Location of the GRD file'
        write(14,530)scale
        530 format(a20)
        write(14,*)'_______ Scale for visualization'
        write(14,*)	amp_min,amp_max
        write(14,*)'_______ scale diapason'

        write(14,59)
        write(14,*)	3
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,560)ver
        560 format('..\..\..\TMP_files\vert\ztr_',a2,'.dat')
        write(14,*)'_______ Location of the DAT file'
        write(14,*)	1
        write(14,*)'_______ Symbol (1: circle, 2: square)'
        write(14,*)	5
        write(14,*)'_______ Size of dots in pixels'
        write(14,*)	250,0,0
        write(14,*)'_______ RGB color'

        write(14,59)
        write(14,*)	3
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,760)ver
        760 format('..\..\..\TMP_files\vert\stat_',a2,'.dat')
        write(14,*)'_______ Location of the DAT file'
        write(14,*)	2
        write(14,*)'_______ Symbol (1: circle, 2: square)'
        write(14,*)	7
        write(14,*)'_______ Size of dots in pixels'
        write(14,*)	0,0,250
        write(14,*)'_______ RGB color'


        write(14,59)
        write(14,*)	2
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,2142)ver
        2142 format('..\..\..\TMP_files\vert\topo_',a2,'.bln')
        write(14,*)'_______ Location of the BLN file'
        write(14,*)	3
        write(14,*)'_______ Thickness of line in pixels'
        write(14,*)	0,0,0
        write(14,*)'_______ RGB color'



        close(14)


        i=system('layers.exe')
        !if(ips.eq.2) stop


541	continue

    end do

    open(14,file='config.txt')
    write(14,*)npix_x,npix_y
    write(14,*)'_______ Size of the picture in pixels (nx,ny)'
    write(14,*)0,dist
    write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
    write(14,*)-zmax,-zmin
    write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
    write(14,*)tick_x,tick_y
    write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
    write(14,251)ar,md,it,ver
    251 format('..\..\..\PICS\',a8,'\',a8,'\IT',a1,'\ver_vpvs',a2,'.png')
    write(14,*)'_______ Path of the output picture'
    write(14,3531)	iver, iver
    3531 format(' Vp/Vs ratio, section=',i1,'A - ',i1,'B')
    write(14,*)'_______ Title of the plot on the upper axe'
    write(14,*)	1
    write(14,*)'_______ Number of layers'

    259 format('********************************************')

    write(14,259)
    write(14,*)	1
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,252)it,ver
    252 format('..\..\..\TMP_files\vert\vpvs_',a1,a2,'.grd')
    write(14,*)'_______ Location of the GRD file'
    write(14,253)scale_vpvs
    253 format(a20)
    write(14,*)'_______ Scale for visualization'
    write(14,*)	vpvs_min,vpvs_max
    write(14,*)'_______ scale diapason'

    write(14,259)
    write(14,*)	3
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,256)ver
    256 format('..\..\..\TMP_files\vert\ztr_',a2,'.dat')
    write(14,*)'_______ Location of the DAT file'
    write(14,*)	1
    write(14,*)'_______ Symbol (1: circle, 2: square)'
    write(14,*)	4
    write(14,*)'_______ Size of dots in pixels'
    write(14,*)	250,0,0
    write(14,*)'_______ RGB color'

    write(14,59)
    write(14,*)	3
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,276)ver
    276 format('..\..\..\TMP_files\vert\stat_',a2,'.dat')
    write(14,*)'_______ Location of the DAT file'
    write(14,*)	2
    write(14,*)'_______ Symbol (1: circle, 2: square)'
    write(14,*)	7
    write(14,*)'_______ Size of dots in pixels'
    write(14,*)	0,0,250
    write(14,*)'_______ RGB color'

    write(14,59)
    write(14,*)	2
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,2141)ver
    2141 format('..\..\..\TMP_files\vert\topo_',a2,'.bln')
    write(14,*)'_______ Location of the BLN file'
    write(14,*)	3
    write(14,*)'_______ Thickness of line in pixels'
    write(14,*)	0,0,0
    write(14,*)'_______ RGB color'



    close(14)

    i=system('layers.exe')

    deallocate(dvan,aprio_dv,vvv,vtmp,vabs,vabsp)

end do
close(16)

open(2,file='../../../DATA/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
close(2)

open(1,file='../../../DATA/'//ar//'/config.txt')
read(1,*) 
read(1,*) npix_h_x0,npix_h_y0
read(1,*) tick_x,tick_y
close(1)

if(npix_h_x0.ne.0) then
    npix_x=npix_h_x0
else
    npix_x=int(npix_h_y0*((fmap2-fmap1)/(tmap2-tmap1)))
    npix_y=npix_h_y0
end if

if(npix_h_y0.ne.0) then
    npix_y=npix_h_y0
else
    npix_y=int(npix_h_x0*((tmap2-tmap1)/(fmap2-fmap1)))
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
write(14,*)tick_x,tick_y
write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
write(14,151)ar,md,it,ps,lv
151 format('..\..\..\PICS\',a8,'\',a8,'\IT',a1,'\hor_profiles',a1,a2,'.png')
write(14,*)'_______ Path of the output picture'
if(iiips.eq.1) then
	write(14,5331)	izzz
5331 format(' P anomalies, depth=',i4,' km')
else
	write(14,5332)	izzz
5332 format(' S anomalies, depth=',i4,' km')
end if
write(14,*)'_______ Title of the plot on the upper axe'
write(14,*)	1
write(14,*)'_______ Number of layers'

write(14,59)
write(14,*)	1
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,552)ps,it,lv
552 format('..\..\..\TMP_files\hor\dv',a1,a1,a2,'.grd')
write(14,*)'_______ Location of the GRD file'
write(14,555)scale_line
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