USE DFPORT

character*8 ar
real fst(1000),tst(1000),zst(1000),tobkr(500)
integer kodst2(1000),istkr(500),ipskr(500)


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0


open(1,file='data_in/SET.DAT')
read(1,*)fi0,tet0,rad_km
close(1)

! Read the coordinates of the stations
open(1,file='data_in/stat_ft.dat')
open(11,file='data_out/stat_ft.dat')
nst=0
33	read(1,*,end=44)fi,tet,zstat
	write(11,*)fi,tet,zstat
	nst=nst+1
	fst(nst)=fi
	tst(nst)=tet
	zst(nst)=zstat
	goto 33
44	close(1)
close(1)
write(*,*)' Total number of stations in the list:',nst
goto 3

222 write(*,*)' file stat_ft.dat does not exist'
stop

3 continue


nztold=0
nztnew=0

open(1,file='data_in/rays.dat')
open(11,file='data_out/rays.dat')

22	continue
    read(1,*,end=21)fzt,tzt,zzt,nkrat
    nztold=nztold+1

    call sfdec(fzt,tzt,0,xzt,yzt,z,fi0,tet0)
    rad=sqrt(xzt**2 + yzt**2)

    if(rad.le.rad_km) write(11,*)fzt,tzt,zzt,nkrat
    if(rad.le.rad_km) nztnew=nztnew+1
    do i=1,nkrat
        nline=nline+1
        read(1,*)ips,ist,tobs
        if(rad.le.rad_km) write(11,*)ips,ist,tobs
    end do
    goto 22

21  close(1)
close(11)
write(*,*)' nztold=',nztold,' nztnew=',nztnew

stop




open(11,file='../../../DATA/'//ar//'/inidata/sources_ini.dat')


nline=0
nray=0
nsrc=0
nst2=0
fztav=0
tztav=0
fmin=9999999
fmax=-9999999
tmin=9999999
tmax=-9999999
np=0
ns=0
! Read the sources:
992	continue
    read(1,*,end=991)fini,tini,zold,nkrat
    write(11,*)fini,tini,zold
    nline=nline+1

    nsrc=nsrc+1
    fztav=fztav+fini
    tztav=tztav+tini

    if(fini.lt.fmin) fmin=fini
    if(fini.gt.fmax) fmax=fini
    if(tini.lt.tmin) tmin=tini
    if(tini.gt.tmax) tmax=tini

    do i=1,nkrat
        nline=nline+1
        read(1,*,end=224,err=224)ips,ist,tobs
        if(ips.ne.1.and.ips.ne.2) goto 225
        if(ist.le.0.or.ist.gt.nst) goto 225
        if(ips.eq.1) np=np+1
        if(ips.eq.2) ns=ns+1

        nray=nray+1
        if(nst2.ne.0) then
            do ist2=1,nst2
                if(ist.eq.kodst2(ist2)) goto 5
            end do
        end if
        nst2=nst2+1
        kodst2(nst2)=ist

5       continue
    end do
    goto 992

991 close(1)
close(11)


write(*,*)' File rays.dat contains lines: ',nline
write(*,*)' Number of events:',nsrc
write(*,*)' Number of picks:',nray
write(*,*)' including P-data:',np,' and S-data:',ns
write(*,*)' Number of stations involved:',nst2

write(31,*)' File rays.dat contains lines: ',nline
write(31,*)' Number of events:',nsrc
write(31,*)' Number of picks:',nray
write(31,*)' including P-data:',np,' and S-data:',ns
write(31,*)' Number of stations involved:',nst2

faver=0
taver=0
open(11,file='../../../DATA/'//ar//'/inidata/stat_actual.dat')

do ist2=1,nst2
    ist=kodst2(ist2)
    fstat2=fst(ist)
    tstat=tst(ist)
    zstat=zst(ist)
    faver=faver+fstat2
    taver=taver+tstat
    write(11,*)fstat2,tstat,zstat
    if(fstat2.lt.fmin) fmin=fstat2
    if(fstat2.gt.fmax) fmax=fstat2
    if(tstat.lt.tmin) tmin=tstat
    if(tstat.gt.tmax) tmax=tstat
end do
close(11)
faver=faver/nst2
taver=taver/nst2

write(*,*)' center of the station network: fi=',faver,' tet=',taver
write(*,*)' fmin=',fmin,' fmax=',fmax
write(*,*)' tmin=',tmin,' tmax=',tmax
write(31,*)' center of the station network: fi=',faver,' tet=',taver
write(31,*)' fmin=',fmin,' fmax=',fmax
write(31,*)' tmin=',tmin,' tmax=',tmax

cost=cos(taver*per)
npix_y=700
npix_x=((fmax-fmin)/(tmax-tmin)) * cost * npix_y


fztav=fztav/nsrc
tztav=tztav/nsrc

write(*,*)' average source point: fi=',fztav,' tet=',tztav
write(31,*)' average source point: fi=',fztav,' tet=',tztav
close(31)

goto 229


stop

224 write(*,*)' Problem in file rays.dat around line:',nline
stop
225 write(*,*)' Problem in file rays.dat around line:',nline
write(*,*)' ips=',ips,' ist=',ist
stop


229 continue

i=system('mkdir ..\..\..\PICS\'//ar)
i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')


open(14,file='config.txt')
write(14,*)npix_x,npix_y
write(14,*)'_______ Size of the picture in pixels (nx,ny)'
write(14,*)fmin,fmax
write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
write(14,*)tmin,tmax
write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
write(14,*)0.25,0.25
write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
write(14,151)ar
151 format('..\..\..\PICS\',a8,'\ini_data.png')
write(14,*)'_______ picture file'
write(14,611)
611 format(' Sources and stations in the initial catalogue')
write(14,*)'_______ Title of the plot on the upper axe'
write(14,*)	1
write(14,*)'_______ Number of layers'

59 format('********************************************')

write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,51)ar
51 format('..\..\..\DATA\',a8,'\inidata\sources_ini.dat')
write(14,*)'_______ Location of the DAT file'
write(14,*)	1
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	4
write(14,*)'_______ Size of dots in pixels'
write(14,*)	255,0,0
write(14,*)'_______ RGB color'


write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,52)ar
52 format('..\..\..\DATA\',a8,'\inidata\stat_actual.dat')
write(14,*)'_______ Location of the DAT file'
write(14,*)	2
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	5
write(14,*)'_______ Size of dots in pixels'
write(14,*)	0,0,255
write(14,*)'_______ RGB color'

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


stop
end