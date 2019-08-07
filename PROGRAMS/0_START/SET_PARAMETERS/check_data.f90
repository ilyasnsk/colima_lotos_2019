USE DFPORT

character*4 dsaa/'DSAA'/
character*8 ar,md
character*100 line100
real fstat2(1000),tstat(1000),zstat(1000)
integer kodst2(1000)
real histo(1000),steploc(4),resloc(4)

allocatable plotray(:,:)

common/pi/pi,per


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0
rz=6371

md='MODEL_01'

open(1,file='SET.SET')
read(1,*)nstep
read(1,*)pmarge
read(1,*)nsmap
read(1,*)nspar
read(1,*)percmin
read(1,*)nkrmin
close(1)


open(1,file='../../../model.dat',status='old',err=1221)
read(1,'(a8)',err=1221)ar		! code of the area
close(1)
goto 2
1221 write(*,*)' file model.dat does not exist'
stop
2 continue
write(*,*)' area=',ar

i=system('mkdir ..\..\..\DATA\'//ar//'\MODEL_01')



open(31,file='../../../DATA/'//ar//'/inidata/data_info.txt')

! Read the coordinates of the stations
open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat',status='old',err=1222)
nst=0
33	read(1,*,end=44)fst,tst,zst
	nst=nst+1
	fstat2(nst)=fst
	tstat(nst)=tst
	zstat(nst)=zst
	goto 33
44	close(1)
close(1)
write(*,*)' Total number of stations in the list:',nst
write(31,*)' Total number of stations in the list:',nst
goto 3

1222 write(*,*)' file stat_ft.dat does not exist'
stop

3 continue


open(1,file='../../../DATA/'//ar//'/inidata/rays.dat',status='old',err=1223)
goto 4

1223 write(*,*)' file rays.dat does not exist'
stop

4 continue

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
1992	continue
    read(1,*,end=1991)fini,tini,zold,nkrat
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
        read(1,*,end=1224,err=1224)ips,ist,tobs
        if(ips.ne.1.and.ips.ne.2) goto 1225
        if(ist.le.0.or.ist.gt.nst) goto 1225
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
    goto 1992

1991 close(1)
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
    fst=fstat2(ist)
    tst=tstat(ist)
    zst=zstat(ist)
    faver=faver+fst
    taver=taver+tst
    write(11,*)fst,tst,zst
    if(fst.lt.fmin) fmin=fst
    if(fst.gt.fmax) fmax=fst
    if(tst.lt.tmin) tmin=tst
    if(tst.gt.tmax) tmax=tst
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

goto 1229


1224 write(*,*)' Problem in file rays.dat around line:',nline
stop
1225 write(*,*)' Problem in file rays.dat around line:',nline
write(*,*)' ips=',ips,' ist=',ist
stop


1229 continue

!*************************************************************************
!*************************************************************************
! COMPUTING RAY DENSITY
!*************************************************************************
!*************************************************************************

dfpl=(fmax-fmin)/nstep
dtpl=(tmax-tmin)/nstep
write(*,*)' dfpl=',dfpl,' dtpl=',dtpl
ds0=(dfpl+dtpl)/20.
allocate(plotray(nstep,nstep))

open(1,file='../../../DATA/'//ar//'/inidata/rays.dat',status='old',err=1223)
nray=0
nsrc=0
av_zzt=0
histo=0
nkrmax=0
ray_max=0
zzt_max=0
ray_aver=0
! Read the sources:
482	continue
    read(1,*,end=481)fzt,tzt,zzt,nkrat
    histo(nkrat)=histo(nkrat)+1
    if(nkrat.gt.nkrmax) nkrmax=nkrat
    if(zzt.gt.zzt_max) zzt_max=zzt
    av_zzt=av_zzt+zzt
    nsrc=nsrc+1
    do i=1,nkrat
        nline=nline+1
        read(1,*)ips,ist,tobs
        fst=fstat2(ist); tst=tstat(ist)
        nray=nray+1
        dist=sqrt((fzt-fst)**2+(tzt-tst)**2)
        dist_deg=epic_dist(fzt,tzt, fst,tst)
        dist_km=dist_deg*per*rz
        if(dist_km.gt.ray_max) ray_max=dist_km
        ray_aver=ray_aver+dist_km
        ns=dist/ds0
        ds=dist/ns
        dfdd=(fst-fzt)/dist
        dtdd=(tst-tzt)/dist
        !if(nray.lt.26283) cycle
        !write(*,*)nray
        do is=1,ns
            sss=(is-0.5)*ds
            fff=fzt+dfdd*sss
            ttt=tzt+dtdd*sss
            do iff=1,nstep
                f1=fmin+(iff-1)*dfpl
                f2=fmin+iff*dfpl
                if((fff-f1)*(fff-f2).le.0.) exit
            end do
            do itt=1,nstep
                t1=tmin+(itt-1)*dtpl
                t2=tmin+itt*dtpl
                if((ttt-t1)*(ttt-t2).le.0.) exit
            end do
            if(iff.le.0.or.iff.gt.nstep) cycle
            if(itt.le.0.or.itt.gt.nstep) cycle
            !write(*,*)' fff=',fff,iff,f1,f2
            !write(*,*)' ttt=',ttt,itt,t1,t2
            plotray(iff,itt)=plotray(iff,itt)+ds
        end do
    end do
    goto 482
481 close(1)
close(11)
av_zzt=av_zzt/nsrc
ray_aver=ray_aver/nray
write(*,*)' MAX depth=',zzt_max,' av_zzt=',av_zzt
write(*,*)' MAX distance=',ray_max,' average distance=',ray_aver

pmax=0
ptot=0
do iff=1,nstep
    do itt=1,nstep
        ptot=ptot+plotray(iff,itt)
        if(plotray(iff,itt).gt.pmax)pmax=plotray(iff,itt)
    end do
end do

histot=0
open(1,file='../../../DATA/'//ar//'/inidata/histo.bln')
write(1,*)nkrmax
do i=1,nkrmax
    a=i
    histot=histot+histo(i)*i
    write(1,*)a,histo(i)
end do
close(1)
write(*,*)' histot=',histot,' nkrmax=',nkrmax

hcur=0
do i=1,nkrmax
    hcur=hcur+i*histo(i)
    if(i.lt.nkrmin) cycle
    write(*,*)' i=',i,' hcur=',hcur,histot*percmin/100
    if(hcur.gt.histot*percmin/100) exit
end do
nkr_lim=i
write(*,*)' nkr_lim=',nkr_lim

!write(*,*)' pmax=',pmax,' ptot=',ptot


!***************************************************************
! Determining limits of the study area
!***************************************************************

porog=ptot*pmarge/100.

pcur=0
do iff=1,nstep
    dp=0
    do itt=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    fmin1=fmin+(iff-1)*dfpl
    exit
end do
!write(*,*)' fmin=',fmin,' fmin1=',fmin1

pcur=0
do iff=nstep,1,-1
    dp=0
    do itt=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    fmax1=fmin+iff*dfpl
    exit
end do
!write(*,*)' fmax=',fmax,' fmax1=',fmax1

pcur=0
do itt=1,nstep
    dp=0
    do iff=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    tmin1=tmin+(itt-1)*dtpl
    exit
end do
!write(*,*)' tmin=',tmin,' tmin1=',tmin1

pcur=0
do itt=nstep,1,-1
    dp=0
    do iff=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    tmax1=tmin+itt*dtpl
    exit
end do

open(14,file='../../../TMP_files/rays/area.bln')
write(14,*)5
write(14,*)fmin1,tmin1
write(14,*)fmin1,tmax1
write(14,*)fmax1,tmax1
write(14,*)fmax1,tmin1
write(14,*)fmin1,tmin1
close(14)

!write(*,*)' tmax=',tmax,' tmax1=',tmax1
fi0=(fmax1+fmin1)/2; tet0=(tmax1+tmin1)/2;
!write(*,*)' pmax=',pmax,' ptot=',ptot

!***************************************************************
! Determining limits of WELL RESOLVED area
!***************************************************************
porog=ptot*10/100.

pcur=0
do iff=1,nstep
    dp=0
    do itt=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    fmin11=fmin+(iff-1)*dfpl
    exit
end do
!write(*,*)' fmin=',fmin,' fmin1=',fmin1

pcur=0
do iff=nstep,1,-1
    dp=0
    do itt=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    fmax11=fmin+iff*dfpl
    exit
end do
!write(*,*)' fmax=',fmax,' fmax1=',fmax1

pcur=0
do itt=1,nstep
    dp=0
    do iff=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    tmin11=tmin+(itt-1)*dtpl
    exit
end do
!write(*,*)' tmin=',tmin,' tmin1=',tmin1

pcur=0
do itt=nstep,1,-1
    dp=0
    do iff=1,nstep
        dp=dp+plotray(iff,itt)
    end do
    pcur=pcur+dp
    if(pcur.lt.porog) cycle
    tmax11=tmin+itt*dtpl
    exit
end do
!write(*,*)' tmax=',tmax,' tmax1=',tmax1

open(14,file='../../../TMP_files/rays/area2.bln')
write(14,*)5
write(14,*)fmin11,tmin11
write(14,*)fmin11,tmax11
write(14,*)fmax11,tmax11
write(14,*)fmax11,tmin11
write(14,*)fmin11,tmin11
close(14)

!***********************************************
! PARAMETERS FOR sethor.dat
!***********************************************
open(11,file='../../../DATA/'//ar//'/sethor.dat')
nmaps=1
write(11,201)nmaps
201 format(i2,'  Number of depth sections for visualization (integer)')
zmap=10
if((tmax1-tmin1).lt.2.)zmap=5
if((tmax1-tmin1).lt.0.5)zmap=1
write(11,202)zmap
202 format(f5.1,'  Depths of sections (real)')
fmap1=fmin1-fi0; fmap2=fmax1-fi0; dfmap=(fmax1-fmin1)/nsmap
tmap1=tmin1-tet0; tmap2=tmax1-tet0; dtmap=(tmax1-tmin1)/nsmap
write(11,203)fmap1,fmap2,dfmap,tmap1,tmap2,dtmap
203 format(2f6.2,f7.4,2f6.2,f7.4,'  fmap1,fmap2,dfmap,tmap1,tmap2,dtmap (real)')
ydist=rz*(tmax1-tmin1)*per
dspar=ydist/nspar
write(11,204)dspar*2
204 format(f4.0,'  Distance to the nearest node (real)')
write(11,205)0
205 format(i1,'  Smoothing (integer)')
write(11,206)1
206 format(i1,'  Visulalization of events (1-yes, 0-no)')
close(1)

!***********************************************
! PARAMETERS FOR setver.dat
!***********************************************
open(11,file='../../../DATA/'//ar//'/setver.dat')
nvert=2
write(11,211)nvert
211 format(i2,'  Number of vertical sections for visualization (integer)')
write(11,212)fmin11,tmin11,fmax11,tmax11
write(11,212)fmin11,tmax11,fmax11,tmin11
212 format(4f10.5,'  Ends of the profile (lon,lat, degrees)')
write(11,213)dspar*2
213 format(f4.0,'  Distance from profile for events,km (real)')
call SFDEC(fmin11,tmin11,0,x1,y1,Z,fi0,tet0)
call SFDEC(fmax11,tmax11,0,x2,y2,Z,fi0,tet0)
dist_km=sqrt((x2-x1)**2+(y2-y1)**2)
write(*,*)' dist_km=',dist_km
ds=dist_km/200
write(11,214)ds
214 format(f4.1,'  Step for visualizing the profile, km (real)')
zver=dist_km/1.5
if(dist_km.lt.300) smark=50
if(dist_km.lt.100) smark=20
if(dist_km.lt.50) smark=10
if(dist_km.lt.20) smark=2
if(dist_km.lt.10) smark=1
!if(av_zzt.lt.zver) zver=av_zzt
zup=-3
write(11,215)zup,zver,ds
215 format(3f7.1,'  Depth limits and step dz, km (real)')
write(11,2151)smark
2151 format(f5.0,'  Marks for indication of position of section')
write(11,216)dspar*2
216 format(f4.0,'  Distance to the nearest node (real)')
write(11,217)0
217 format(i1,'  Smoothing (integer)')
write(11,218)1
218 format(i1,'  Visulalization of events (1-yes, 0-no)')
close(11)

!***********************************************
! PARAMETERS FOR config.txt
!***********************************************
open(11,file='../../../DATA/'//ar//'/config.txt')
write(11,*)'******* MAP VIEW *********'
npix_x=600
costet=cos(tet0*per)
apix_y=npix_x * ((tmax1-tmin1)/(fmax1-fmin1)) / costet
npix_y=apix_y
write(11,221)npix_x,npix_y
221 format(i5,i5,' size in pixels for horizontal section (int)')
tick=(fmax1-fmin1)/5
write(*,*)' tick=',tick
if((fmax1-fmin1).gt.5.)then
    write(11,2221)tick,tick
    2221 format(2f4.0,' ticks on axes for horizontal section, deg (real)')
else if((fmax1-fmin1).gt.2.)then
    write(11,2222)tick,tick
    2222 format(2f5.1,' ticks on axes for horizontal section, deg (real)')
else
    write(11,2223)tick,tick
    2223 format(2f6.2,' ticks on axes for horizontal section, deg (real)')
end if
write(11,*)'******* VERTICAL SECTION *********'
npix_x=0
npix_y=600
write(11,223)npix_x,npix_y
223 format(i5,i5,' size in pixels for vertical section (int)')
stick=dist_km/5
write(11,224)stick,stick
224 format(2f4.0,' size in pixels for vertical section (int)')
write(11,*)'******* PLOTS WITH 1D VELOCITIES *********'
write(11,*)'500 500			size in pixels for the 1D models'
write(11,*)'2 9				Limits of P and S velocity distribution'
write(11,*)'-60 5			Depth limit'
write(11,*)'0.5 10			ticks on axes for 1D velocity plot'
write(11,*)'******* SCALES *********'
write(11,*)'blue_red.scl	scale for velocity anomalies'
write(11,*)'-10 10			diapason for velocity anomalies, %'
write(11,*)'blue_brown.scl	scale for Vp/Vs'
write(11,*)'1.6 1.9			diapason for Vp/Vs'
write(11,*)'rainbow_small.scl	scale for absolute Vp'
write(11,*)'4.0 8.5			diapason for absolute Vp'
write(11,*)'rainbow_small.scl	scale for absolute Vp'
write(11,*)'2.5 5	'
close(11)

! Limits 

!***********************************************
! PARAMETERS FOR MAJOR_PARAM.DAT
!***********************************************
open(11,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
write(11,*)'********************************************************'
line100='GENERAL INFORMATION :'
write(11,'(a100)')line100
write(11,*)'1	KEY 1: REAL; KEY 2: SYNTHETIC ' 
write(11,*)'1	KEY 1: Vp and Vs; KEY 2: Vp and Vp/Vs  '
write(11,*)'0	KEY 0: all data, KEY 1: odd events, KEY 2: even events'
write(11,*)'0	Ref. model optimization (0-no; 1-yes)'
write(11,*)''
write(11,*)''
write(11,*)'********************************************************'
line100='AREA_CENTER :'
write(11,'(a100)')line100
if((fmax1-fmin1).gt.5.)then
    write(11,2311)fi0,tet0
    2311 format(2f6.0,' Center of conversion to XY, degrees (real)')
else if((fmax1-fmin1).gt.2.)then
    write(11,2312)fi0,tet0
    2312 format(2f7.1,' Center of conversion to XY, degrees (real)')
else 
    write(11,2313)fi0,tet0
    2313 format(2f8.2,' Center of conversion to XY, degrees (real)')
end if
write(11,*)
write(11,*)
write(11,*)'********************************************************'
line100='ORIENTATIONS OF GRIDS :'
write(11,'(a100)')line100
write(11,*)'4				number of grids'
write(11,*)'0 22 45 67		orientations'
write(11,*)
write(11,*)
write(11,*)'********************************************************'
line100='1D LOCATION KEY :'
write(11,'(a100)')line100
lokey=1
if(dist_km.lt.50) lokey=2
if(dist_km.lt.200.and.av_zzt.gt.30) lokey=2
write(11,232)lokey
232 format(i1,'    1: using reference table (large areas); (int)')
write(11,*)'     2: using straight lines (small areas with high relief)'
write(11,*)
write(11,*)
write(11,*)'********************************************************'
line100='INVERSION PARAMETERS :'
write(11,'(a100)')line100
write(11,*)'100 		LSQR iterations, iter_max'
write(11,*)'1 1.		Weights for P and S models '
write(11,*)'0.6 1.3		HORIZONTAL smoothing (P, S)'
write(11,*)'0.6 1.3		VERTICAL smoothing (P, S)'
write(11,*)'0.0 0.0		AMPLITUDE damping (P, S)'
write(11,*)''
write(11,*)''
write(11,*)'1.0001   1.0001	weight of the station corrections (P and S)'
write(11,*)'5.0	wzt_hor'
write(11,*)'5.0	wzt_ver'
write(11,*)'2.0	wzt_time'
write(11,*)
write(11,*)
write(11,*)'********************************************************'
write(11,*)'Parameters for location in 1D model using reference table'
write(11,*)'and data selection:'
write(11,*)'********************************************************'
line100='LIN_LOC_PARAM :'
write(11,'(a100)')line100
write(11,233)nkr_lim
233 format(i2,'		Minimal number of records')
write(11,*)'100		km, maximum distance to nearest station'
write(11,*)'1.5		S max resid with respect to P max resid '
write(11,*)'100		dist_limit=100	: within this distance the weight is equal'
write(11,*)'1		n_pwr_dist=1	: power for decreasing of W with distance'
write(11,*)'30		ncyc_av=10	'
write(11,*)''
write(11,*)'! For output:'
write(11,*)'30		bad_max=30		: maximal number of outliers'
if(dist_km.gt.200) then
    dtdd=0.1; dist_dt=15
else if(dist_km.gt.50) then
    dtdd=0.1; dist_dt=10
else
    dtdd=0.2; dist_dt=2.5
end if
write(*,*)' dist_km=',dist_km,' dtdd=',dtdd,' dist_dt=',dist_dt
write(11,234)dtdd
234 format(f5.2,' dt_dd, maximal dt/distance for threshold of residuals')
write(11,235)dist_dt
235 format(f5.2,' dist_dt, distance limit (dtmax=dt_dd*dist_dt)')
write(11,*)''
if(nsrc.gt.10000) then
    nprint=500
else if(nsrc.gt.4000) then
    nprint=200
else if(nsrc.gt.1000) then
    nprint=100
else 
    nprint=50
end if
write(*,*)' nsrc=',nsrc,' nprint=',nprint
write(11,236)nprint
236 format(i5,' Frequency for output printing')

if(dist_km.gt.300) then
    nit_loc=3
    steploc(1)=20;resloc(1)=8
    steploc(2)=5;resloc(2)=4
    steploc(3)=1;resloc(3)=2
else if(dist_km.gt.150) then
    nit_loc=3
    steploc(1)=10;resloc(1)=5
    steploc(2)=3;resloc(2)=3
    steploc(3)=1;resloc(3)=2
else if(dist_km.gt.50) then
    nit_loc=3
    steploc(1)=5;resloc(1)=4
    steploc(2)=2;resloc(2)=2
    steploc(3)=0.5;resloc(3)=1
else
    nit_loc=2
    steploc(1)=2;resloc(1)=2
    steploc(2)=0.3;resloc(2)=1
end if
write(11,*)
write(11,237)nit_loc
237 format(i1,' Number of different grids for searching events')
do i=1,nit_loc
    write(11,*)'_______________________________________________________'
    write(11,238)steploc(i),steploc(i),steploc(i)
    238 format(3f6.1,' dx,dy,dz, spacing of grid for searching (km)')
    write(11,239)0.
    239 format(f5.1,' res_loc1, lower limit for location (for LT residuals, W=1)')
    write(11,240)resloc(i)
    240 format(f5.1,' res_loc1, upper limit for location (for GT residuals, W=0)')
    write(11,241)2.
    241 format(f5.1,' w_P_S_diff (+ causes better coherency of P and S)')
end do
write(11,*)
write(11,*)
write(11,*)'********************************************************'
write(11,*)'Parameters for location in 3D model using bending tracing'
write(11,*)'********************************************************'
line100='LOC_PARAMETERS:'
write(11,'(a100)')line100
write(11,*)'! Parameters for BENDING:'
ds_integr=ray_max/400
if(ray_aver/ds_integr.gt.100) ds_integr=ray_aver/100
write(*,*)' Integration step for ray tracing: ds_integr=',ds_integr
write(11,242)ds_integr
242 format(f6.2,' Integration step for ray tracing')
segm_min=ray_aver/5
write(11,243)segm_min
243 format(f6.2,' Min segment for bending')
bend_min=ray_aver/100
write(11,244)bend_min
244 format(f6.2,' Min value of bending')
bend_max=ray_aver/3
write(11,245)bend_max
245 format(f6.2,' Max value of bending')
write(11,*)''
write(11,*)'! Parameters for location'
write(11,*)'30		dist_limit=100	: within this distance the weight is equal'
write(11,*)'1		n_pwr_dist=1	: power for decreasing of W with distance'
write(11,*)'30		ncyc_av=10	'
write(11,*)''
write(11,*)'0.		res_loc1=0.2	: lower limit for location (for LT residuals, W=1)'
write(11,246)resloc(nit_loc)
246 format(f4.1,' res_loc1, upper limit for location (for GT residuals, W=0)')
write(11,*)'2.		w_P_S_diff=2 (+ causes better coherency of P and S)'
write(11,247)steploc(1)
247 format(f4.1,' Max step for shifting event after inversion')
write(11,248)steploc(nit_loc)
248 format(f4.1,' Min step for shifting event')
write(11,*)''
write(11,249)nprint
249 format(i5,' Frequency for output printing')

call SFDEC(fmax1,tmax1,0.,xsize,ysize,Z,fi0,tet0)
size1=sqrt(xsize**2+ysize**2)
call SFDEC(fmax11,tmax11,0.,xsize,ysize,Z,fi0,tet0)
size0=sqrt(xsize**2+ysize**2)
ds_grid0=size0/30.
nds=int(ds_grid0*10); ads=nds
ds_grid=ads/10.
size=ds_grid*30
z_grid=size/1.5
write(*,*)' xsize=',xsize,' ysize=',ysize,' size1=',size1
write(*,*)' size=',size,' ds_grid=',ds_grid
write(11,*)
write(11,*)
write(11,*)'********************************************************'
write(11,*)'Parameters for grid construction'
write(11,*)'********************************************************'
line100='GRID_PARAMETERS:'
write(11,'(a100)')line100

write(11,250)-size,size,ds_grid
250 format(3f6.1,' Limits and spacing for grid construction along X')
write(11,251)-size,size,ds_grid
251 format(3f6.1,' Limits and spacing for grid construction along Y')
write(11,252)-5.,z_grid,ds_grid
252 format(3f6.1,' Limits and spacing for grid construction along Z')
write(11,*)'0.05 100.0	!plotmin, plotmax= maximal ray density, relative to average'
write(11,*)
write(11,*)
write(11,*)'********************************************************'
write(11,*)'Parameters for 3D model with regular grid'
write(11,*)'********************************************************'
line100='3D_MODEL PARAMETERS:'
write(11,'(a100)')line100
write(11,253)-size,size,ds_grid
253 format(3f6.1,' Limits and spacing for regular grid: X')
write(11,254)-size,size,ds_grid
254 format(3f6.1,' Limits and spacing for regular grid: Y')
write(11,255)-5.,z_grid,ds_grid
255 format(3f6.1,' Limits and spacing for regular grid: Z')
write(11,256)ds_grid*10
256 format(f6.1,' distance to the nearest node')
write(11,*)'0		Smoothing factor1'

if(ray_max.gt.200) then
    ddd=500
    if(ray_max.gt.500.) ddd=ray_max*1.1
    zzz=500
    if(zzt_max.gt.500.) zzz=zzt_max*1.1
    dstep=1
else
    ddd=200
    if(ray_max.gt.200.) ddd=ray_max*1.1
    zzz=200
    if(zzt_max.gt.200.) zzz=zzt_max*1.1
    dstep=0.5
end if

write(11,*)'********************************************************'
write(11,*)'Parameters for calculation of the reference table:'
write(11,*)'********************************************************'
line100='REF_PARAM:'
write(11,'(a100)')line100
write(11,*)''
write(11,257)dstep
257 format(f6.1,' step for epicentral distance to compute table')
write(11,258)zzz
258 format(f6.1,' max depth to stop computing table')
write(11,259)ddd
259 format(f6.1,' max distance to stop computing table')
if(ray_max.gt.200) then
    write(11,*)'3		number of depth steps'
    write(11,*)'-5 1	depth, step'
    write(11,*)'40 2	depth, step'
    write(11,*)'200 5	depth, step'
else
    write(11,*)'3		number of depth steps'
    write(11,*)'-5 0.5	depth, step'
    write(11,*)'20 1	depth, step'
    write(11,*)'50 3	depth, step'
end if
write(11,260)zzz
260 format(f6.1,' lower depth limit for events')
close(11)

i=system('copy ..\..\..\COMMON\various\ref_start.dat ..\..\..\DATA\'//ar//'\'//md//'\ref_start.dat')

open(11,file='../../../all_areas.dat')
write(11,*)'1: name of the area (any 8 characters)'
write(11,*)'2: name of the model (any 8 characters)'
write(11,*)'3: number of iterations'
write(11,*)'**********************************************'
write(11,261)ar,md,3
261 format(a8,1x,a8,1x,i1)

stop


i=system('mkdir ..\..\..\TMP_files\rays')

open(14,file='../../../TMP_files/rays/plotray_ini.grd')
write(14,'(a4)')dsaa
write(14,*)nstep,nstep
write(14,*)fmin,fmax
write(14,*)tmin,tmax
write(14,*)0,pmax
do itet=1,nstep
	write(14,*)(plotray(ifi,itet),ifi=1,nstep)
end do
close(14)


stop




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