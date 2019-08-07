character*8 ar,md,line
character*12 file_att
character*1 itt,ps
character*5 stac,stacod(10000),stacur(1000)
real fstat(10000),tstat(10000),zstat(10000)
real attkr(1000),gdkr(1000)
integer istkr(1000)

common/ray_param/ds_ini,ds_part_min,val_bend_min,bend_max0
common/ray/ nodes,xray(1000),yray(1000),zray(1000)
common/center/fi0,tet0
common/pi/pi,per
one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0


open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
close(1)

write(*,*)' RAY TRACING for ATTENUATION: ar=',ar,' md=',md

nzt1=0
open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays0.dat')
746 read(1,*,end=846)xzt,yzt,zzt,nkrat
    nzt1=nzt1+1
    if(nkrat.gt.0) then
        do i=1,nkrat
            read(1,*)ips,ist,tob
        end do
    end if
    goto 746
846 close(1)

nzt2=0
np=0; ns=0

open(3,file='../../../DATA/'//ar//'/'//md//'/data/rays_att.dat')
747 read(3,*,end=847)natt
!write(*,*)natt
    nzt2=nzt2+1
    if(natt.gt.0) then
        do i=1,natt
            read(3,*)ips,ist,tob
            !write(*,*)ips,ist,tob
            if(ips.eq.11) np=np+1
            if(ips.eq.12) ns=ns+1
        end do
    end if
    goto 747

847 close(3)
    
open(11,file='../../../DATA/'//ar//'/'//md//'/data/num_ps_att.dat')
write(11,*)np,ns
close(11)

    
if(nzt1.ne.nzt2) then
    write(*,*)' PROBLEM: number of events for vel and att is not equal: nzt1=',nzt1,' nzt2=',nzt2
    stop
end if
write(*,*)' np=',np,' ns=',ns


!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=573)line
	if(line.eq.'ATTENUAT') goto 574
end do
573 continue
write(*,*)' cannot find ATTENUATION PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
574 read(1,*)npmin,nsmin
read(1,*)iter
close(1)

write(itt,'(i1)')iter

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

call read_vref(ar,md)
call read_3D_mod_v(ar,md,iter)

open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
nst=0
3 read(1,*,end=4)fst,tst,zst,stac
    !write(*,*)fst,tst,zst,stac
    if(stac.eq.'     ') goto 4
    nst=nst+1
    fstat(nst)=fst; tstat(nst)=tst; zstat(nst)=zst; stacod(nst)=stac;

    goto 3
4 close(1)
write(*,*)' nst=',nst

!******************************************************************

do ips12=1,2
    if(ips12.eq.1.and.np.lt.npmin) cycle
    if(ips12.eq.2.and.ns.lt.nsmin) cycle

    write(ps,'(i1)')ips12
    write(*,*)' ps=',ps
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays'//itt//'.dat')
    open(2,file='../../../DATA/'//ar//'/'//md//'/data/rays_att.dat')

    open(11,file='../../../DATA/'//ar//'/'//md//'/data/atten'//ps//'.dat')
    open(12,file='../../../DATA/'//ar//'/'//md//'/data/ray_paths_att'//ps//'.dat',form='binary')
    open(21,file='../../../TMP_files/rays/ray_paths_att'//ps//'.bln')
    open(22,file='../../../TMP_files/rays/srces_att'//ps//'.dat')
    open(23,file='../../../TMP_files/rays/stat_att'//ps//'.dat')
    open(38,file='../../../TMP_files/rays/att_vs_time'//ps//'.dat')
    nray=0
    nzt1=0
    nzt2=0
    att_aver=0

    7   continue
    read(1,*,end=8)xzt,yzt,zzt,nkrat
    !write(*,*)xzt,yzt,zzt,nkrat
    nzt1=nzt1+1
    if(nkrat.gt.0) then
        do i=1,nkrat
            read(1,*)ips,ist,tob
        end do
    end if

    read(2,*)natt
    if(natt.gt.0) then
        do i=1,natt
            read(2,*)ips,ist,tob
            if(ips12.eq.1.and.ips.ne.11) cycle
            if(ips12.eq.2.and.ips.ne.12) cycle
            
            atten=1000./tob

            nray=nray+1
            fst=fstat(ist); tst=tstat(ist)
            call SFDEC(fst,tst,0.,xst,yst,Z,fi0,tet0)
            call trace_bending(xzt,yzt,zzt,xst,yst,zst,ips12,	tout)

            write(11,*)xzt,yzt,zzt,xst,yst,zst,atten,tout
            call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
            call decsf(xst,yst,0.,fi0,tet0,fst2,tst2,h)

            write(38,*)tout,atten

            att_aver=att_aver+(atten/tout)

            if(mod(nray,50).eq.0) write(*,*)nray,' ist=',ist,' atten/tout=',atten/tout,' tout=',tout

            write(21,*)2
            write(21,*)fzt,tzt
            write(21,*)fst,tst

            write(22,*)fzt,tzt,zzt
            write(23,*)fst,tst,zst

            write(12)nodes
            !write(*,*)' nodes=',nodes
            !write(*,*)' ztr=',xzt,yzt,zzt
            do inod=1,nodes
                write(12)xray(inod),yray(inod),zray(inod)
                !write(*,*)xray(inod),yray(inod),zray(inod)
            end do
            !write(*,*)' sta=',xst,yst,zst
            !stop
        end do
    end if

    goto 7

    8 close(1);close(2)
    close(11);close(12);close(21);close(22);close(23);close(38)

    att_aver=att_aver/nray

    write(*,*)' ips=',ips12,' nray=',nray,' att_aver=',att_aver

end do



stop
end