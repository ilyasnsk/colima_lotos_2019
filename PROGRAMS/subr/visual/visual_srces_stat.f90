subroutine visual_srces_stat(ar,md,it)
character*1 it
character*8 ar,md,line
common/center/fi0,tet0

!******************************************************************
open(1,file='../../data/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=553)line
	if(line.eq.'AREA_CEN') goto 554
end do
553 continue
write(*,*)' cannot find AREA CENTER in MAJOR_PARAM.DAT!!!'
pause
554 read(1,*)fi0,tet0
close(1)


open(1,file='../../data/'//ar//'/config.txt')
read(1,*) 
read(1,*) npix_x,npix_y
read(1,*) tick_x,tick_y
close(1)

open(2,file='../../data/'//ar//'/sethor.dat')
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
write(14,151)ar,md,it
151 format('..\..\PICS\',a8,'\',a8,'\sources',a1,'.png')
write(14,*)'_______ Path of the output picture'
if(it.eq.'a') then
	write(14,*)	' Stations (blue) and initial sources (red)'
else 
	write(14,*)	' Stations (blue) and sources (red) after iteration ',it
end if
write(14,*)'_______ Title of the plot on the upper axe'
write(14,*)	1
write(14,*)'_______ Number of layers'

59 format('********************************************')

write(14,59)
write(14,*)	2
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,55)ar
55 format('..\..\DATA\',a8,'\map\coastal_line.bln')
write(14,*)'_______ Location of the BLN file'
write(14,*)	3
write(14,*)'_______ Thickness of line in pixels'
write(14,*)	0,0,0
write(14,*)'_______ RGB color'

write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,156)ar
156 format('..\..\DATA\',a8,'\inidata\stat_ft.dat')
write(14,*)'_______ Location of the DAT file'
write(14,*)	2
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	8
write(14,*)'_______ Size of dots in pixels'
write(14,*)	0,0,255
write(14,*)'_______ RGB color'


write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
if(it.eq.'a') then
	write(14,158)ar
	158 format('..\..\DATA\',a8,'\inidata\srces.dat')
else
	write(14,157)ar,md,it
	157 format('..\..\DATA\',a8,'\',a8,'\data\srces',a1,'.dat')
end if
write(14,*)'_______ Location of the DAT file'
write(14,*)	1
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	5
write(14,*)'_______ Size of dots in pixels'
write(14,*)	255,0,0
write(14,*)'_______ RGB color'



close(14)

i=system('copy ..\..\COMMON\visual_exe\visual.exe layers.exe')

i=system('layers.exe')


return
end