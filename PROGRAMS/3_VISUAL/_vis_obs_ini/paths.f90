! This program visualize the rays corresponding to current calculations
! (from 'tmp' directory)
! Output:
! FIG_files/rays/rays_ver.dat : projection of rays in vertical section (as points)
! FIG_files/rays/ztr_ver.dat : projection of event coordinates in vertical section
! FIG_files/rays/rays_hor.dat : rays in a defined layer, z_up - z_low (as points)
! FIG_files/rays/ztr_hor.dat : sources in map view

USE DFPORT

character*8 ar,md,line
common/pi/pi,per

one=1
pi=asin(one)*2.e0
per=pi/180.e0
rz=6371.

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
close(1)
write(*,*)' ar=',ar,' md=',md

!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=553)line
	if(line.eq.'AREA_CEN') goto 554
end do
553 continue
write(*,*)' cannot find AREA CENTER in MAJOR_PARAM.DAT!!!'
pause
554 read(1,*)fi0,tet0
write(*,*)fi0,tet0
close(1)


! Read the coordinates of the stations
open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
nst=0
33	read(1,*,end=44)fi,tet,zstat
    nst=nst+1
    goto 33
44	close(12)
close(1)

open(1,file='../../../DATA/'//ar//'/inidata/rays.dat')
open(11,file='../../../DATA/'//ar//'/inidata/events_ini.dat')
izt=0
nray=0
nrp=0
nrs=0

! Read the sources:
992	continue
    read(1,*,end=991)fzt,tzt,zzt,nkrat
    write(11,*)fzt,tzt,zzt,nkrat
    izt=izt+1
    do i=1,nkrat
        read(1,*)ips,ist,tobs
        nray=nray+1
        if(ips.eq.1)nrp=nrp+1
        if(ips.eq.2)nrs=nrs+1
    end do
goto 992
991 close(1)
close(11)

write(*,*)' Number of stations=',nst
write(*,*)' Number of events:',izt
write(*,*)' Total rays:',nray
write(*,*)' P-rays:',nrp
write(*,*)' S-rays:',nrs

open(21,file='../../../DATA/'//ar//'/inidata/INFO.TXT')
write(21,*)' Number of stations=',nst
write(21,*)' Number of events:',izt
write(21,*)' Total rays:',nray
write(21,*)' P-rays:',nrp
write(21,*)' S-rays:',nrs
close(21)


!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************
!*****************************************************************

open(1,file='../../../DATA/'//ar//'/config.txt')
read(1,*) 
read(1,*) npix_x,npix_y
read(1,*) tick_x,tick_y
close(1)

open(2,file='../../../DATA/'//ar//'/sethor.dat')
read(2,*)  
read(2,*)   
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
close(2)


open(14,file='config.txt')
write(14,*)npix_x,npix_y
write(14,*)'_______ Size of the picture in pixels (nx,ny)'
write(14,*)fmap1+fi0,fmap2+fi0
write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
write(14,*)tmap1+tet0,tmap2+tet0
write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
write(14,*)tick_x,tick_y
write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
write(14,151)ar,md
151 format('..\..\..\PICS\',a8,'\',a8,'\obs_system_ini_hor.png')
write(14,*)'_______ Path of the output picture'
write(14,*)	' Stations (blue) and initial sources (red)'
write(14,*)'_______ Title of the plot on the upper axe'
write(14,*)	1
write(14,*)'_______ Number of layers'

59 format('********************************************')

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

write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,158)ar
158 format('..\..\..\DATA\',a8,'\inidata\events_ini.dat')
write(14,*)'_______ Location of the DAT file'
write(14,*)	1
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	5
write(14,*)'_______ Size of dots in pixels'
write(14,*)	255,0,0
write(14,*)'_______ RGB color'


write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,156)ar
156 format('..\..\..\DATA\',a8,'\inidata\stat_ft.dat')
write(14,*)'_______ Location of the DAT file'
write(14,*)	2
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	8
write(14,*)'_______ Size of dots in pixels'
write(14,*)	0,0,255
write(14,*)'_______ RGB color'


close(14)

i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')

i=system('layers.exe')


stop
end