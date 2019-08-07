USE DFPORT

character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*1 ps
character*2 lv
real hlev(100)
character*20 scale_dv,scale_vpvs
integer line1_rgb(3),line2_rgb(3),kdot_rgb(3)

allocatable v_ini(:,:),v_abs(:,:,:)
common/center/fi0,tet0
common/pi/pi,per
common/keys/key_ft1_xy2

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

open(1,file='../../../model.dat')
read(1,'(a8)')ar	! synthetic model
read(1,'(a8)')md	! synthetic model
close(1)

i=system('mkdir ..\..\..\TMP_files\hor')


write(*,*)' Synthetic model in horizontal sections'
write(*,*)' ar=',ar,' md=',md

key_preview=0
open(1,file='../../../preview_key.txt')
read(1,*,end=771)key_preview
771 close(1)

if(key_preview.ne.0) then
        i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\SYN')
    open(1,file='../../../DATA/'//ar//'/config.txt')
    read(1,*) 
    read(1,*) npix_x0,npix_y0
    read(1,*)tick_x,tick_y
    read(1,*) 
    read(1,*) 
    read(1,*) 
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
    close(1)
    i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')
    i=system('copy ..\..\..\COMMON\scales_scl\'//scale_dv//' '//scale_dv)
    i=system('copy ..\..\..\COMMON\scales_scl\'//scale_vpvs//' '//scale_vpvs)


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
read(1,*)
read(1,*)
read(1,*,end=441,err=441)key_ft1_xy2
441 close(1)
!******************************************************************

!write(*,*)' key_ft1_xy2=',key_ft1_xy2
!pause


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

open(2,file='../../../data/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
close(2)

if(npix_x0.ne.0) then
    npix_x=npix_x0
else
    npix_x=int(npix_y0*((fmap2-fmap1)/(tmap2-tmap1)))
    npix_y=npix_y0
end if

if(npix_y0.ne.0) then
    npix_y=npix_y0
else
    npix_y=int(npix_x0*((tmap2-tmap1)/(fmap2-fmap1)))
    npix_x=npix_x0
end if

write(*,*)' npix_x=',npix_x,' npix_y=',npix_y


dfmap=dfmap/2.
dtmap=dtmap/2.

nfmap=int_best((fmap2-fmap1)/dfmap+1.)
ntmap=int_best((tmap2-tmap1)/dtmap+1.)
write(*,*)' nfmap=',nfmap,' ntmap=',ntmap


allocate(v_ini(nfmap,ntmap),v_abs(nfmap,ntmap,2))

DO ilev=1,nlev
    zzz=hlev(ilev)
    write(lv,'(i2)')ilev
    write(*,*)' ilev=',ilev,' zzz=',zzz
    v_abs=0
    do ips=1,2
        write(ps,'(i1)')ips
        vref=vrefmod(zzz,ips)

        v_ini=0

        do itet=1,ntmap
            ttt=(itet-1)*dtmap+tmap1+tet0
            !write(*,*)' itet=',itet,' ttt=',ttt
            !ttt=-7.35
            do ifi=1,nfmap
                fff=(ifi-1)*dfmap+fmap1+fi0
                relief=h_lim(fff,ttt)

                !fff=-60; ttt=13

                if (zzz.lt.relief) then
                    v_ini(ifi,itet)=-999
                    cycle
                end if
                if(key_ft1_xy2.eq.1) then
                    call SFDEC(fff,ttt,0.,xxx,yyy,Z,fi0,tet0)
                else
                    xxx=fff
                    yyy=ttt
                end if
                dv=anomaly(xxx,yyy,zzz,ips)
                !write(*,*)' fi=',fff,' tet=',ttt,' dv=',dv
                !write(*,*)' x=',xxx,' y=',yyy,' z=',zzz,' dv=',dv
                !pause

                if (zzz.lt.h_lim(fff,ttt)) dv=0
                v_ini(ifi,itet)=dv
                v_abs(ifi,itet,ips)=vref*(1+dv/100)


                !pause
                !pause
            end do
        end do
  
       open(14,file='../../../TMP_files/hor/syn'//ps//'_'//lv//'.xyz')
       do itet=1,ntmap
            ttt=(itet-1)*dtmap+tmap1+tet0
            do ifi=1,nfmap
                fff=(ifi-1)*dfmap+fmap1+fi0
                write(14,*)fff,ttt,v_ini(ifi,itet)
            end do
        end do
        close(14)


        open(14,file='../../../TMP_files/hor/syn'//ps//'_'//lv//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nfmap,ntmap
        write(14,*)fmap1+fi0,fmap2+fi0
        write(14,*)tmap1+tet0,tmap2+tet0
        write(14,*)-999,999
        do itet=1,ntmap
            write(14,*)(v_ini(ifi,itet),ifi=1,nfmap)
        end do
        close(14)

        if(key_preview.eq.0) cycle


        !*********************************************************
        !*********************************************************
        !*********************************************************
        !*********************************************************
        !*********************************************************
        !*********************************************************

        open(14,file='config.txt')
        write(14,*)npix_x,npix_y
        write(14,*)'_______ Size of the picture in pixels (nx,ny)'
        write(14,*)fmap1+fi0,fmap2+fi0
        write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
        write(14,*)tmap1+tet0,tmap2+tet0
        write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
        write(14,*)tick_x,tick_y
        write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
        write(14,51)ar,md,ps,lv
        51 format('..\..\..\PICS\',a8,'\',a8,'\SYN\hor_syn',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        izzz=zzz
        if(ips.eq.1) then
        write(14,3331)	izzz
        3331 format(' P anomalies, depth=',i5)
        else
        write(14,3332)	izzz
        3332 format(' S anomalies, depth=',i5)
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'

        59 format('********************************************')

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,52)ps,lv
        52 format('..\..\..\TMP_files\hor\syn',a1,'_',a2,'.grd')
        write(14,*)'_______ Location of the GRD file'
        write(14,53)scale_dv
        53 format(a20)
        write(14,*)'_______ Scale for visualization'
        write(14,*)	dv_min,dv_max
        write(14,*)'_______ scale diapason'


        open(1,file='../../../data/'//ar//'/map/polit_bound.bln',status='old',err=491)
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

        open(1,file='../../../data/'//ar//'/map/coastal_line.bln',status='old',err=492)
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

        close(14)

        i=system('layers.exe')


    end do


! IMAGE THE VP/VS RATIO:


    open(14,file='../../../TMP_files/hor/syn_vpvs_'//lv//'.xyz')
    do itet=1,ntmap
        ttt=(itet-1)*dtmap+tmap1+tet0
        do ifi=1,nfmap
            fff=(ifi-1)*dfmap+fmap1+fi0
            write(14,*)fff,ttt,v_abs(ifi,itet,1)/v_abs(ifi,itet,2)
        end do
    end do
    close(14)



    open(14,file='../../../TMP_files/hor/syn_vpvs_'//lv//'.grd')
    write(14,'(a4)')dsaa
    write(14,*)nfmap,ntmap
    write(14,*)fmap1+fi0,fmap2+fi0
    write(14,*)tmap1+tet0,tmap2+tet0
    write(14,*)-999,999
    do itet=1,ntmap
        write(14,*)(v_abs(ifi,itet,1)/v_abs(ifi,itet,2),ifi=1,nfmap)
    end do
    close(14)

    if(key_preview.eq.0) cycle


    !*********************************************************
    !*********************************************************
    !*********************************************************
    !*********************************************************
    !*********************************************************
    !*********************************************************

    open(14,file='config.txt')
    write(14,*)npix_x,npix_y
    write(14,*)'_______ Size of the picture in pixels (nx,ny)'
    write(14,*)fmap1+fi0,fmap2+fi0
    write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
    write(14,*)tmap1+tet0,tmap2+tet0
    write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
    write(14,*)tick_x,tick_y
    write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
    write(14,451)ar,md,lv
    451 format('..\..\..\PICS\',a8,'\',a8,'\SYN\hor_syn_vpvs',a2,'.png')
    write(14,*)'_______ Path of the output picture'
    izzz=zzz
    write(14,3431)	izzz
    3431 format(' Vp/Vs ratio, depth=',i5)
    write(14,*)'_______ Title of the plot on the upper axe'
    write(14,*)	1
    write(14,*)'_______ Number of layers'

    write(14,59)
    write(14,*)	1
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,452)lv
    452 format('..\..\..\TMP_files\hor\syn_vpvs_',a2,'.grd')
    write(14,*)'_______ Location of the GRD file'
    write(14,453)scale_vpvs
    453 format(a20)
    write(14,*)'_______ Scale for visualization'
    write(14,*)	vpvs_min,vpvs_max
    write(14,*)'_______ scale diapason'


    open(1,file='../../../data/'//ar//'/map/polit_bound.bln',status='old',err=691)
    close(1)
    write(14,59)
    write(14,*)	2
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,454)ar
    454 format('..\..\..\DATA\',a8,'\map\polit_bound.bln')
    write(14,*)'_______ Location of the BLN file'
    write(14,*)	3
    write(14,*)'_______ Thickness of line in pixels'
    write(14,*)	100,0,0
    write(14,*)'_______ RGB color'
    691 continue

    open(1,file='../../../data/'//ar//'/map/coastal_line.bln',status='old',err=692)
    close(1)
    write(14,59)
    write(14,*)	2
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,455)ar
    455 format('..\..\..\DATA\',a8,'\map\coastal_line.bln')
    write(14,*)'_______ Location of the BLN file'
    write(14,*)	3
    write(14,*)'_______ Thickness of line in pixels'
    write(14,*)	0,0,0
    write(14,*)'_______ RGB color'
    692 continue

    close(14)

    i=system('layers.exe')


end do

stop
end