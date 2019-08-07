USE DFPORT

character*8 ar,md,line,md_ini,ar_ini
character*1 ps,itt,it0,rm,gr,it_ini
common/ray/ nodes,xray(2000),yray(2000),zray(2000)
real xstn(1000),ystn(1000),zstn(1000),statinv(2,1000)
common/nanom/n_anomaly
common/pi/pi,per
common/center/fi0,tet0
common/ray_param/ds_ini,ds_part_min,val_bend_min,bend_max0
common/noise_2/kod_noise,ar_ini,md_ini,it_ini,red_ps(2)
common/keys/key_ft1_xy2


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
close(1)

ar_ini=ar


i=system('mkdir ..\..\..\TMP_files\1D_mod')
i=system('mkdir ..\..\..\DATA\'//ar//'\'//md//'\data')

write(*,*)' Forward modeling: ar=',ar,' md=',md

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

!******************************************************************
kod_ini1_mod2=1
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=253)line
	if(line.eq.'SYNTHETI') goto 254
end do
253 continue
write(*,*)' cannot find SYNTHETIC MODELING PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
254 read(1,'(a8)',err=257)md_ini
if(md_ini.ne.'inidata ') kod_ini1_mod2=2
read(1,*,err=257)iter_ini
write(it_ini,'(i1)')iter_ini
257 close(1)
!******************************************************************

write(*,*)' ar_ini=',ar_ini,' md_ini=',md_ini,' it_ini=',it_ini

call prepare_noise(ar,md)

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
close(1)

!******************************************************************

open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=583)line
	if(line.eq.'LOC_PARA') goto 584
end do
583 continue
write(*,*)' cannot find LOC_PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
584 continue

read(1,*)	! Parameters for bending:
read(1,*)ds_ini
read(1,*)ds_part_min
read(1,*)val_bend_min
read(1,*)bend_max0
close(1)
w_qual=1
!******************************************************************


call read_topo(ar)
call read_vref(ar,md)
call read_anom(ar,md)

!goto 1313
!		zzt=60
!		xzt=-130
!		yzt=20
!		xst=140
!		yst=20
!		zst=0
!		ips=1
!		open(21,file='../../../TMP_files/rays/ray_work.bln')
!		do xst=-90.,160.,5.
!			call trace_bending(xzt,yzt,zzt,xst,yst,zst,ips,	tout)
!			write(*,*)xst,' tout=',tout
!			write(21,*)nodes
!			do i=1,nodes
!				x=xray(i)
!				y=yray(i)
!				z=zray(i)
!				dv = dv_board(x,y,z,ips)
!				v = velocity(x,y,z,ips)
!				!write(*,*)x,y,z,dv
!				write(21,*)x+1.5*per*6371.,-z,dv
!			end do
!		end do
!		close(21)
!		stop
!
!1313 continue


! Read the coordinates of the stations
open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
nst=0
33	read(1,*,end=44)fi,tet,zstat
	call SFDEC(fi,tet,0.,X,Y,Z,fi0,tet0)
	nst=nst+1
	xstn(nst)=x
	ystn(nst)=y
	zstn(nst)=zstat
	goto 33
44	continue
close(1)
!write(*,*)' nst=',nst

if(kod_ini1_mod2.eq.2) then
    open(1,file='../../../DATA/'//ar_ini//'/'//md_ini//'/data/rays'//it_ini//'.dat',status='old',err=571)
    goto 572
571 write(*,*)' Cannot find file ../../../DATA/',ar_ini,'/',md_ini,'/data/rays',it_ini,'.dat'
    stop
    572 continue
else 
    open(1,file='../../../DATA/'//ar//'/inidata/rays.dat')
end if
open(11,file='../../../DATA/'//ar//'/'//md//'/data/rays_syn.dat')

nzt=0
nray=0
nr=0
disp_tot1=0
disp_tot2=0
disp=0
dd=0
dp=0
ds=0
np=0
ns=0
21	continue

    if(kod_ini1_mod2.eq.2) then
        read(1,*,end=22)xzt,yzt,zzt,nkrat
        call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
    else
        read(1,*,end=22)fzt,tzt,zzt,nkrat
        call SFDEC(fzt,tzt,0.,xzt,yzt,Z,fi0,tet0)
    end if
    nzt=nzt+1
    !if(nzt.gt.200) goto 22
    !	call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
    write(11,*)fzt,tzt,zzt,nkrat

    do ikrat=1,nkrat
    if(kod_ini1_mod2.eq.2) then
        read(1,*)ips,ist,tobs,tref
        resid=tobs-tref
    else
        read(1,*)ips,ist,tobs
    end if
    xst=xstn(ist)
    yst=ystn(ist)
    zst=zstn(ist)
    nray=nray+1
    call trace_bending(xzt,yzt,zzt,xst,yst,zst,ips,	tout)
    !if(nray.le.867) cycle

    if(kod_noise.eq.1) then
        dt_rand=our_noise(nray,0.0,ips)
    else if(kod_noise.eq.2) then
        dt_rand=resid*red_ps(ips)
    else
        dt_rand=0.
    end if

    !write(*,*)' dt_rand=',dt_rand
    tout=tout+dt_rand
    dd=dd+abs(dt_rand)

    write(11,*)ips,ist,tout
    if(ips.eq.1)np=np+1
    if(ips.eq.2)ns=ns+1


    end do
    ddcur=dd/nray
    if(mod(nzt,10).eq.0)write(*,*)' nzt=',nzt,' np=',np,' ns=',ns,' er=',ddcur
    goto 21
22 close(1)
close(11)

write(*,*)' nzt=',nzt,' nray=',nray,' np=',np,' ns=',ns




stop
end