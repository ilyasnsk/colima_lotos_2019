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

i=system('mkdir ..\..\..\TMP_files\att')


write(*,*)' Synthetic model in horizontal sections'
write(*,*)' ar=',ar,' md=',md

key_preview=1
key_ft1_xy2=1
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
!******************************************************************
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

call read_topo(ar)
call read_anom(ar,md)

open(2,file='../../../data/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
close(2)

cost=cos(tet0*per)
if(npix_x0.ne.0) then
    npix_x=npix_x0
else
    npix_x=int(npix_y0*((fmap2-fmap1)/(tmap2-tmap1))*cost)
    npix_y=npix_y0
end if

if(npix_h_y0.ne.0) then
    npix_y=npix_y0
else
    npix_y=int(npix_x0*((tmap2-tmap1)/(fmap2-fmap1))/cost)
    npix_x=npix_x0
end if

write(*,*)' npix_x=',npix_x,' npix_y=',npix_y


dfmap=dfmap/2.
dtmap=dtmap/2.

nfmap=int_best((fmap2-fmap1)/dfmap+1.)
ntmap=int_best((tmap2-tmap1)/dtmap+1.)
write(*,*)' nfmap=',nfmap,' ntmap=',ntmap


allocate(v_ini(nfmap,ntmap))
do ips=1,2
    write(ps,'(i1)')ips
    DO ilev=1,nlev
        zzz=hlev(ilev)
        write(lv,'(i2)')ilev
        write(*,*)' ilev=',ilev,' zzz=',zzz

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
                dv=anomaly(xxx,yyy,zzz,ips+2)
                !write(*,*)' fi=',fff,' tet=',ttt,' dv=',dv
                !write(*,*)' x=',xxx,' y=',yyy,' z=',zzz,' dv=',dv
                !pause

                if (zzz.lt.h_lim(fff,ttt)) dv=0
                v_ini(ifi,itet)=dv


                !pause
                !pause
            end do
        end do
  

        open(14,file='../../../TMP_files/att/syn_hor'//ps//lv//'.grd')
        write(14,'(a4)')dsaa
        write(14,*)nfmap,ntmap
        write(14,*)fmap1+fi0,fmap2+fi0
        write(14,*)tmap1+tet0,tmap2+tet0
        write(14,*)-999,999
        do itet=1,ntmap
            write(14,*)(v_ini(ifi,itet),ifi=1,nfmap)
        end do
        close(14)

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
        51 format('..\..\..\PICS\',a8,'\',a8,'\SYN\att_hor',a1,a2,'.png')
        write(14,*)'_______ Path of the output picture'
        izzz=zzz
        if(ips.eq.1) then
        write(14,3331)	izzz
        3331 format(' P attenuation, depth=',i5)
        else
        write(14,3332)	izzz
        3332 format(' S attenuation, depth=',i5)
        end if
        write(14,*)'_______ Title of the plot on the upper axe'
        write(14,*)	1
        write(14,*)'_______ Number of layers'

        59 format('********************************************')

        write(14,59)
        write(14,*)	1
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,52)ps,lv
        52 format('..\..\..\TMP_files\att\syn_hor',a1,a2,'.grd')
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
end do

stop
end