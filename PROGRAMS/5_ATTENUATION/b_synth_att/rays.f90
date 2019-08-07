USE DFPORT

character*8 ar,md,line,md_ini,ar_ini
character*1 ps,itt,it0,rm,gr,it_ini
common/ray/ nodes,xray(2000),yray(2000),zray(2000)
real xstn(1000),ystn(1000),zstn(1000),statinv(2,1000)
common/nanom/n_anomaly
common/pi/pi,per
common/center/fi0,tet0
common/ray_param/ds_ini,ds_part_min,val_bend_min,bend_max0
common/keys/key_ft1_xy2


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)ips		! code of the area
close(1)

ar_ini=ar
write(ps,'(i1)') ips


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
if(ips.eq.2) att_aver=val21; add_att=val22
close(1)
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
goto 258
257 write(*,*)' the model name does not exist in MAJOR_PARAM'
stop
258 close(1)
!******************************************************************
write(*,*)' ar=',ar,' md_ini=',md_ini,' md=',md,' ps=',ps   
i=system('copy ..\..\..\DATA\'//ar//'\'//md_ini//'\data\ray_paths_att'//ps//&
'.dat ..\..\..\DATA\'//ar//'\'//md//'\data\ray_paths_att'//ps//'.dat')

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


call read_anom(ar,md)



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
write(*,*)' nst=',nst

open(1,file='../../../DATA/'//ar//'/'//md_ini//'/data/atten'//ps//'.dat')
open(2,file='../../../DATA/'//ar//'/'//md_ini//'/data/ray_paths_att'//ps//'.dat',form='binary')

open(11,file='../../../DATA/'//ar//'/'//md//'/data/atten'//ps//'.dat')

nzt=0
nray=0

aver=0
21	continue

    read(1,*,end=22)xzt,yzt,zzt,xst,yst,zst,atten,ttt
    read(2)nodes
    if(nodes.eq.0)goto 21
    do i=1,nodes
        read(2)xray(i),yray(i),zray(i)
        !write(*,*)xray(i),yray(i),zray(i)
    end do
    nray=nray+1
    dasyn=0.
    sss=0
    do ipt=2,nodes
        x1=xray(ipt-1)
        y1=yray(ipt-1)
        z1=zray(ipt-1)
        x2=xray(ipt)
        y2=yray(ipt)
        z2=zray(ipt)
        ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
        xmid=(x1+x2)/2.
        ymid=(y1+y2)/2.
        zmid=(z1+z2)/2.
        an=anomaly(xmid,ymid,zmid,ips+2)
        !an=5
        !write(*,*)xmid,ymid,zmid,an
        dasyn=dasyn+an*ds
        sss=sss+ds
    end do
    dasyn=dasyn/sss
    a_ref = add_att + att_aver*ttt
    a_rand=our_noise(nray,0.0,1)
    !write(*,*)' dt_rand=',dt_rand
    a_syn = a_ref*(1 + dasyn/100.) + a_rand
    aver=aver+abs(dasyn)
    if(mod(nray,10).eq.0) write(*,*)nray,' a_ref, dasyn, a_syn=',a_ref,dasyn,a_syn
    !pause

    !att_abs=att_aver+asyn

    write(11,*)xzt,yzt,zzt,xst,yst,zst,a_syn,ttt

    goto 21
22 close(1)
close(11)

avcur=aver/nray
write(*,*)' nray=',nray,' aver=',avcur




stop
end