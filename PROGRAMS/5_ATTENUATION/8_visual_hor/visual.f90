USE DFPORT
! HORIZONTAL !!!
! ATTENUATION !!!!!!!

character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*2 lv
character*1 ps

real hlev(100)

character*20 scale_line, scale_line2
allocatable dvan(:,:),vvv(:,:)

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
read(1,*)ips
close(1)
write(ps,'(i1)')ips

write(*,*)' ar=',ar,' md=',md,' ps=',ps

i=system('mkdir ..\..\..\TMP_files\att')
i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\ATT')

key_preview=1
open(1,file='../../../DATA/'//ar//'/config.txt')
read(1,*) 
read(1,*)npix_x0,npix_y0
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
read(1,*)scale_line
read(1,*)a1,a2,amp_min,amp_max
read(1,*)scale_line2
read(1,*)vpvs_min,vpvs_max
close(1)

i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')
i=system('copy ..\..\..\COMMON\scales_scl\'//scale_line//' '//scale_line)
i=system('copy ..\..\..\COMMON\scales_scl\'//scale_line2//' '//scale_line2)




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

open(2,file='../../../DATA/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
read(2,*) smaxx
close(2)

cost=cos(tet0*per)
if(npix_x0.ne.0) then
    npix_x=npix_x0
else
    npix_x=int(npix_y0*((fmap2-fmap1)/(tmap2-tmap1))*cost)
    npix_y=npix_y0
end if

if(npix_y0.ne.0) then
    npix_y=npix_y0
else
    npix_y=int(npix_x0*((tmap2-tmap1)/(fmap2-fmap1))/cost)
    npix_x=npix_x0
end if

write(*,*)' npix_x=',npix_x,' npix_y=',npix_y

call read_topo(ar)

nfmap=int_best((fmap2-fmap1)/dfmap+1.)
ntmap=int_best((tmap2-tmap1)/dtmap+1.)
write(*,*)' nfmap=',nfmap,' ntmap=',ntmap
allocate(dvan(nfmap,ntmap),vvv(nfmap,ntmap))

DO ilev=1,nlev
    zzz=hlev(ilev)
    write(lv,'(i2)')ilev
    write(*,*)' ilev=',ilev,' zzz=',zzz
    dvan=0
    vvv=0
    do igr=ngr1,ngr2
        call prepare_model_att(ar,md,igr,ips)

        do itet=1,ntmap
            ttt=(itet-1)*dtmap+tmap1+tet0
            !ttt=tet0
            do ifi=1,nfmap
                fff=(ifi-1)*dfmap+fmap1+fi0
                !fff=fi0
                !fff=110.87
                dv=0
                www=0
                call att_1_grid(fff,ttt,zzz,smaxx, att,www)
                !write(*,*)' att=',att,' www=',www
                !stop
                zlim_up=h_lim(fff,ttt)
                if (zzz.lt.zlim_up) www=0

                dvan(ifi,itet)=dvan(ifi,itet)+att*www
                vvv(ifi,itet)=vvv(ifi,itet)+www
            end do
        end do
    end do
    amax=-99999
    amin=999999
    do ifi=1,nfmap
        do itet=1,ntmap
            attcur=-999.
            if (vvv(ifi,itet).gt.w_limit) then
	            att=dvan(ifi,itet)/vvv(ifi,itet)
                attcur=100.*att/att_aver
                if(attcur.gt.amax) amax=attcur
                if(attcur.lt.amin) amin=attcur
            end if
            dvan(ifi,itet)=attcur
        end do
    end do
       
    open(14,file='../../../TMP_files/att/att_hor'//ps//lv//'.grd')
    write(14,'(a4)')dsaa
    write(14,*)nfmap,ntmap
    write(14,*)fmap1+fi0,fmap2+fi0
    write(14,*)tmap1+tet0,tmap2+tet0
    write(14,*)amin,amax
    do itet=1,ntmap
	    write(14,*)(dvan(ifi,itet),ifi=1,nfmap)
    end do
    close(14)

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
	51 format('..\..\..\PICS\',a8,'\',a8,'\ATT\hor_att',a1,a2,'.png')
	write(14,*)'_______ Path of the output picture'
	izzz=zzz
	if(ips.eq.1) then
		write(14,3331)	izzz
	3331 format(' P attenuation, depth=',i3,' km')
	else
		write(14,3332)	izzz
	3332 format(' S attenuations, depth=',i3,' km')
	end if
	write(14,*)'_______ Title of the plot on the upper axe'
	write(14,*)	1
	write(14,*)'_______ Number of layers'

	59 format('********************************************')

	write(14,59)
	write(14,*)	1
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,52)ps,lv
	52 format('..\..\..\TMP_files\att\att_hor',a1,a2,'.grd')
	write(14,*)'_______ Location of the GRD file'
	write(14,53)scale_line
	53 format(a20)
	write(14,*)'_______ Scale for visualization'
	write(14,*)	amp_min,amp_max
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

	close(14)


	i=system('layers.exe')

end do

stop
end