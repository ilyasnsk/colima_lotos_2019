USE DFPORT

character*8 ar,md,line
character*1 ps,itt,it0,gr
PARAMETER (nkrmax=500)
real xstn(5000),ystn(5000),zstn(5000),statinv(2,1000)

real dt3(500), w_qual(10)
real aaa(500,4),atmp(500,4),bbb(500),btmp(500),xxx(4)
real www(4),vvv(4,4),dttest(500)
real tob_best(500),trf_best(500),dtold(500)

allocatable rays_all(:,:,:),rays_best(:,:,:),nod_all(:),nod_best(:)

real xvert1(20),yvert1(20),xvert2(20),yvert2(20)
character*2 ver
real xmark(200,20),ymark(200,20),smark(200,20)
integer nmark(100)
character*5 stacod(500),stac


common/refmod/nrefmod,zref(600),vref(600,2)
common/pi/pi,per
common/ray_param/ds_ini,ds_part_min,val_bend_min,bend_max0

common/loc_param/wgs,res_loc1,res_loc2,dist_limit,n_pwr_dist,ncyc_av,w_P_S_diff
common/krat/nkrat,istkr(500),tobkr(500),ipskr(500),qualkr(500),trfkr(500),ngood(500),alkr(500),diskr(500)

common/ray/ nodes,xray(1000),yray(1000),zray(1000)
common/nkr_max/nnkrmax
common/center/fi0,tet0
common/keys/key_ft1_xy2
common/flat_sf/key_flat1

rz=6371
nnkrmax=nkrmax
ifreq=10


one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

r_hor=0.1
r_ver=0.01
r_time=0.1

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)iter		! code of the grid
close(1)

write(itt,'(i1)')iter
write(it0,'(i1)')iter-1

write(*,*)' SOURCE LOCATION: ar=',ar,' md=',md,' it=',itt

!******************************************************************
key_ft1_xy2=1
key_true1=0
key_flat1=1
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
read(1,*,end=441,err=441)key_ft1_xy2    ! 1: geogr, 2: Cartesian, 
read(1,*,end=441,err=441)key_true1      ! 1: true locations exist in separate line
read(1,*,end=441,err=441)key_flat1      ! 1: calculations in flat model, 2: spherical velocity
441 close(1)
!******************************************************************

write(*,*)' key_flat1=',key_flat1


i=system('mkdir ..\..\..\TMP_files\tmp')
i=system('mkdir ..\..\..\TMP_files\loc')

if(key_ft1_xy2.eq.1) then
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
else
    fi0=0; tet0=0
end if


call read_3D_mod_v(ar,md,iter-1)
call read_z_lim(ar,md)
call read_ini_model(ar,md)
call read_topo(ar)



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
write(*,*)nornt

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
read(1,*)
read(1,*)	! Parameters for location:
read(1,*)dist_limit	!=100
read(1,*)n_pwr_dist	!=1
read(1,*)ncyc_av	!=10
read(1,*)
read(1,*)res_loc1
read(1,*)res_loc2
read(1,*)w_P_S_diff
read(1,*)stepmax
read(1,*)stepmin
read(1,*)
read(1,*)ifreq
close(1)
w_qual=1

call read_vref(ar,md)

key_info=0
open(3,file='../../../DATA/'//ar//'/inidata/event_info.dat',status='old',err=46)
key_info=1
close(3)
46  continue


! Read the coordinates of the stations
open(2,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
if(key_info.eq.1) open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
i=0
3 i=i+1
    read(2,*,end=4)xstn(i),ystn(i),zstn0
    if(key_info.eq.1)read(1,*)fi,tet,zst,stacod(i)

    zstn(i) = flat_sph(key_flat1,xstn(i),ystn(i),zstn0)

    !write(*,*)xstn(i),ystn(i),zstn0,zstn(i)
    goto 3
4 close(2)
if(key_info.eq.1) close(1)
nst=i-1
!write(*,*)' nst=',nst

statinv=0
if(iter.ne.1) then
    do igr=1,nornt

	write(gr,'(i1)')igr
	open(2,file='../../../DATA/'//ar//'/'//md//'/data/stcor_p_'//it0//gr//'.dat')
	do ist=1,nst
	    read(2,*,end=332)stcor
	    statinv(1,ist) = statinv(1,ist)+stcor
	end do
	close(2)

332    open(2,file='../../../DATA/'//ar//'/'//md//'/data/stcor_s_'//it0//gr//'.dat')
	do ist=1,nst
	    read(2,*,end=333)stcor
	    statinv(2,ist) = statinv(2,ist)+stcor
	end do
	close(2)
    end do
    statinv = statinv / nornt
    333	continue
end if


if(iter.ne.1) then
	do igr=1,nornt
		write(gr,'(i1)')igr
		iun=30+igr
		open(iun,file='../../../DATA/'//ar//'/'//md//'/data/ztcor_'//it0//gr//'.dat')
	end do
end if

!write(*,*)' key_info=',key_info

!goto 1313

open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays'//it0//'.dat')

open(11,file='../../../DATA/'//ar//'/'//md//'/data/rays'//itt//'.dat')
open(41,file='../../../DATA/'//ar//'/'//md//'/data/resid'//itt//'.dat')
open(14,file='../../../DATA/'//ar//'/'//md//'/data/srces'//itt//'.dat')
open(12,file='../../../TMP_files/tmp/ray_paths_'//itt//'.dat',form='binary')
open(21,file='../../../TMP_files/rays/dt_dist_P_'//itt//'.dat')
open(22,file='../../../TMP_files/rays/dt_dist_S_'//itt//'.dat')

if(key_info.eq.1)then
    open(3,file='../../../DATA/'//ar//'/'//md//'/data/info_'//it0//'.dat')
    open(15,file='../../../DATA/'//ar//'/'//md//'/data/info_'//itt//'.dat')
    open(16,file='../../../DATA/'//ar//'/'//md//'/data/rays_full_'//itt//'.dat')
end if


nzt=0
nray=0
dis_tot1=0
dis_tot2=0
ntot=0
err_loc0=0
err_loc=0
nray_s=0
nray_p=0
21	continue
    if(key_true1.eq.1) read(1,*,end=22)xtrue,ytrue,ztrue
    read(1,*,end=22)xini,yini,depth0,nkrat
    !write(*,*)' depth0=',depth0

    if(key_info.eq.1) read(3,*)myr1,mmnt1,mdy1,mhr1,mmn1,sec1

    zini = flat_sph(key_flat1,xini,yini,depth0)

    !write(*,*)' zini=',zini,' depth0=',depth0

    !write(*,*)xini,yini,zini,nkrat
    dx_corr=0
    dy_corr=0
    dz_corr=0
    dt_corr=0

    if(iter.ne.1) then
        do igr=1,nornt
            iun=30+igr
            read(iun,*)dx,dy,dz,dt
            !write(*,*)igr,dx,dy,dz,dt
            dx_corr=dx_corr+dx
            dy_corr=dy_corr+dy
            dz_corr=dz_corr+dz
            dt_corr=dt_corr+dt
        end do
        dx_corr=dx_corr/nornt
        dy_corr=dy_corr/nornt
        dz_corr=dz_corr/nornt
        dt_corr=dt_corr/nornt
       ! write(*,*)' sum:',dx_corr,dy_corr,dz_corr,dt_corr

    end if

    xzt=xini - dx_corr
    yzt=yini - dy_corr
    zzt=zini - dz_corr
    !write(*,*)xzt,yzt,zzt,nkrat


    d_ztr=sqrt(xzt*xzt+yzt*yzt)
    nzt=nzt+1
    !if(zzt.lt.0) zzt=0

    if(key_ft1_xy2.eq.1) then
        call decsf(xzt,yzt,zzt,fi0,tet0,fff,ttt,hhh)
    else
        fff=xxt; ttt=yzt
    end if

    zlim_up=h_lim(fff,ttt)

   ! write(*,*)' zlim_up=',zlim_up,' zzt=',zzt


    if (zzt.le.zlim_up) then
        depth=zlim_up
        zzt = flat_sph(key_flat1,xzt,yzt,depth)
    end if
    !write(*,*)' zlim_up=',zlim_up,' zzt=',zzt

    if(nkrat.eq.0) goto 21
    nkr=0
    do ikrat=1,nkrat
        read(1,*)ips,ist,tobs_old,tref
       ! write(*,*)ist,ips,tobs_old,tref
        nray=nray+1
        tobs=tobs_old - dt_corr - statinv(ips,ist)

        istkr(ikrat)=ist
        ipskr(ikrat)=ips
        tobkr(ikrat)=tobs
        trfkr(ikrat)=tref
        dtold(ikrat)=tobs_old-tref
        qualkr(ikrat)=1
    end do
    

! TEMP !!!!!!!!!!!!!!!!!
!if(nzt.ne.21) goto 21


    !if (d_ztr.gt.d_ztr_max) goto 21

    allocate(rays_all(3,1000,nkrat),rays_best(3,1000,nkrat),nod_all(nkrat),nod_best(nkrat))

    itstep=0
    dstot=0
    dscur=0

    goal_best=0
    step_cur=stepmax
    !write(*,*)' step_cur=',step_cur
    if(mod(nzt,ifreq).eq.0)write(*,*)xzt,yzt,zzt,nkrat

    331	continue
    itstep=itstep+1
    vzt1=velocity(xzt,yzt,zzt,1)
    vzt2=velocity(xzt,yzt,zzt,2)


    !write(*,*)' source:',xzt,yzt,zzt,nkrat
    nk=0
    ddd1=0
    ddd2=0


    do ikrat=1,nkrat
        ist=istkr(ikrat)
        ips=ipskr(ikrat)
        !if(ist.ne.35.or.ips.ne.2) cycle
        xst=xstn(ist)
        yst=ystn(ist)
        zst=zstn(ist)
        dshor=sqrt((xst-xzt)*(xst-xzt)+(yst-yzt)*(yst-yzt))
        dsver=abs(zzt-zst)
        dist=sqrt(dshor*dshor+dsver*dsver)
        diskr(ikrat)=dist

        s0=1/vzt1
        if(ips.ne.1) s0=1/vzt2

        tobs=tobkr(ikrat)
        tref=trfkr(ikrat)
        resid=tobs-tref
        !write(*,*)ikrat, ' ips=',ips,' resid=',resid,' dist=',dist
        !write(*,*)xzt,yzt,zzt

        call trace_bending(xzt,yzt,zzt,xst,yst,zst,ips,	tout)
        !write(*,'(2i6,5f10.3)')ist,nodes,xst,yst,tobs,tout,tobs-tout

        ddd1=ddd1+abs(dtold(ikrat))
        ddd2=ddd2+abs(tobs-tout)
        !write(*,*)' tobs=',tobs,' tout=',tout,trfkr(ikrat),' nodes=',nodes

        nod_all(ikrat)=nodes
        do i=1,nodes
            rays_all(1,i,ikrat)=xray(i)
            rays_all(2,i,ikrat)=yray(i)
            rays_all(3,i,ikrat)=zray(i)
        end do

!if(iter.eq.1) then
!call trace_bending(xzt+2,yzt-3,zzt+4,xst,yst,zst,ips,	tobs)
!tobkr(ikrat)=tobs
!end if

        trfkr(ikrat)=tout
        x1=xray(1)
        y1=yray(1)
        z1=zray(1)
        x2=xray(2)
        y2=yray(2)
        z2=zray(2)


        hor=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
        tot=sqrt(hor*hor+(z2-z1)*(z2-z1))
        !write(*,*)' hor=',hor,' tot=',tot

        cosbe=(x2-x1)/hor
        sinbe=(y2-y1)/hor
        !write(*,*)' dx=',x2-x1,' dy=',y2-y1,' dz=',z2-z1

        cosal=(z2-z1)/tot
        sinal=hor/tot

        px=cosbe*sinal*s0
        py=sinbe*sinal*s0
        pz=cosal*s0
        nk=nk+1
        aaa(nk,1)=px
        aaa(nk,2)=py
        aaa(nk,3)=pz
        aaa(nk,4)=1
        bbb(nk)=tobs-tout



        !write(*,*)' ztr:',xzt,yzt,zzt
        !do i=1,nodes
	        !write(*,*)xray(i),yray(i),zray(i)
        !end do
        !write(*,*)' sta:',xst,yst,zst
        !pause


        nkr=nkr+1

        dt3(ikrat)=tobs-tout
		
        !write(*,*)ips,dist,' dt=',dt3(ikrat)

    end do
    red=100*(ddd1-ddd2)/ddd1
    	!write(*,*)' d1=',ddd1/nkrat,' d2=',ddd2/nkrat,' red=',red


    if(nkr.eq.0) goto 21


    do i1=1,nkrat
        bbb(i1)=dt3(i1)
        if(ipskr(i1).ne.2) cycle
        do i2=1,nkrat
            if(i2.eq.i1) cycle
            if(istkr(i1).ne.istkr(i2)) cycle
            if(ipskr(i2).ne.1) then
	            if(mod(nzt,ifreq).eq.0)write(*,*)' i1=',i1,' ips1=',ipskr(i1),' i2=',i2,' ips2=',ipskr(i2)
            end if
            ipskr(i1)=3

            dt3(i1)=dt3(i1)-dt3(i2)
            bbb(i1)=dt3(i1)
            do ii=1,4
	            aaa(i1,ii)=(aaa(i1,ii)-aaa(i2,ii))
            end do
            exit
        end do
    end do

    call dispers(dt3,	disp3,aver3,nkr3,ank)
!    pause



    do i1=1,nkrat
        !write(*,*)' nkrat=',nkrat,' tobkr(nkrat)=',tobkr(nkrat)
        tobkr(i1)=tobkr(i1) -aver3
        dt3(i1)=tobkr(i1) - trfkr(i1)
        ddd=diskr(i1)
        !write(*,*)' i1=',i1,' ddd=',ddd
        if(ddd.lt.dist_limit) then
	        wdist = 1
        else
	        wdist = (dist_limit/(ddd))**n_pwr_dist
        end if

        ccc=1
        if(ipskr(i1).eq.3) ccc = w_P_S_diff
        if(ipskr(i1).ne.3) bbb(i1) = bbb(i1)-aver3
        do ii=1,4
	        aaa(i1,ii)=aaa(i1,ii)*wdist*ccc
        end do
        bbb(i1) = bbb(i1)*wdist*ccc
    end do

!    do i=1,nkrat
!	    write(*,*)ipskr(i),' t, dt=',trfkr(i),tobkr(i)-trfkr(i)
!    end do
!
    atmp=aaa
    btmp=bbb

    nk=0
    aaa=0
    bbb=0
    do ik=1,nkrat
        if(ipskr(ik).eq.3) ipskr(ik)=2
        if(abs(dt3(ik)).gt.res_loc2) cycle
        nk=nk+1
        bbb(nk)=btmp(ik)
        do i4=1,4
	        aaa(nk,i4)=atmp(ik,i4)
        end do
        !write(*,*)' ist=',istkr(i),' ips=',ipskr(i)
    end do


    if(ank.gt.goal_best) then
        x_best=xzt
        y_best=yzt
        z_best=zzt
        tob_best=tobkr
        trf_best=trfkr
        goal_best=ank
        s_best=dstot
        rays_best=rays_all
        nod_best=nod_all
    else
        dscur=dscur/2.
        dxcur=dxcur/2.
        dycur=dycur/2.
        dzcur=dzcur/2.

        xzt=xzt-dxcur
        yzt=yzt-dycur
        zzt=zzt-dzcur


        if(key_ft1_xy2.eq.1) then
            call decsf(xzt,yzt,zzt,fi0,tet0,fff,ttt,hhh)
        else
            fff=xxt; ttt=yzt
        end if

        zlim_up=h_lim(fff,ttt)
        depth = sph_flat(key_flat1,xzt,yzt,zzt)
        if (depth.le.zlim_up) then
            depth=zlim_up
            zzt = flat_sph(key_flat1,xzt,yzt,depth)
        end if
        !write(*,*)xzt,yzt,zzt


        dstot=dstot-dscur
        !write(*,*)' ***************** dstot=',dstot,' dscur=',dscur
        step_cur=dscur
        if(dscur.gt.stepmin) goto 331

    end if


    do i4=1,4
	    nk=nk+1
	    reg=r_hor
	    if(i4.eq.3)reg=r_ver
	    if(i4.eq.4)reg=r_time
	    aaa(nk,i4)=reg
    end do


    atmp = aaa
    btmp = bbb
    m=nk
    n=4
    mp=nkrmax
    np=4

    call SVDCMP(aaa,M,N,MP,NP,www,vvv)
    call SVBKSB(aaa,www,vvv,M,N,MP,NP,bbb,xxx)
    dt=xxx(4)
    !write(*,*)' dx=',xxx(1),' dy=',xxx(2),' dz=',xxx(3),' dt=',dt
    disp1=0
    disp2=0
    do i=1,nk-4
        dt=0
        do j=1,4
	        dt=dt+atmp(i,j)*xxx(j)
        end do
        disp1=disp1+abs(btmp(i))
        disp2=disp2+abs(btmp(i)-dt)
        !write(*,'(3f8.4)')bbb(i),bbb(i)-dt,dt
    end do
    disp1=disp1/nk
    disp2=disp2/nk
    red=100*(disp1-disp2)/disp1
   !write(*,*)' disp1=',disp1,' disp2=',disp2,' red=',red


    shift0=sqrt(xxx(1)*xxx(1)+xxx(2)*xxx(2)+xxx(3)*xxx(3))
    !write(*,*)' shift0=',shift0

    scale=1.
    if(shift0.gt.step_cur)scale=step_cur/shift0

    dxcur=-xxx(1)*scale
    dycur=-xxx(2)*scale
    dzcur=-xxx(3)*scale

    xnew=xzt+dxcur
    ynew=yzt+dycur
    znew=zzt+dzcur

    !write(*,*)' xnew=',xnew,' ynew=',ynew,' znew=',znew

    if(key_ft1_xy2.eq.1) then
        call decsf(xnew,ynew,0.,fi0,tet0,fff,ttt,hhh)
    else
        fff=xnew; ttt=ynew
    end if
    depth = sph_flat(key_flat1,xnew,ynew,znew)
!    write(*,*)' znew=',znew,' depth=',depth

    limit=0
!    write(*,*)' fff=',fff,' ttt=',ttt
    zlim_up=h_lim(fff,ttt)
!    write(*,*)' zlim_up=',zlim_up
    if (depth.le.zlim_up) then
        depth=zlim_up
        znew = flat_sph(key_flat1,xnew,ynew,depth)
        limit=1
    end if


    deplim = z_lim(fff,ttt)
!    write(*,*)' deplim=',deplim
    if (depth.ge.deplim) then
        depth=deplim
        znew = flat_sph(key_flat1,xnew,ynew,depth)
        limit=1
    end if

!    write(*,*)' znew=',znew
!
!    write(*,*)' limit=',limit
	
    if(limit.eq.1) then
        nk=nk-4
        do i=1,nk
	        aaa(i,3)=1
        end do
        do i4=1,3
            nk=nk+1
            reg=r_hor
            if(i4.eq.3)reg=r_time
            aaa(nk,i4)=reg
        end do
        atmp = aaa
        btmp = bbb
        m=nk
        n=3
        mp=nkrmax
        np=3

        call SVDCMP(aaa,M,N,MP,NP,www,vvv)
        call SVBKSB(aaa,www,vvv,M,N,MP,NP,bbb,xxx)
        dt=xxx(3)
        !write(*,*)' dx=',xxx(1),' dy=',xxx(2),' dt=',dt
        disp1=0
        disp2=0
        do i=1,nk-3
	    dt=0
	    do j=1,3
		    dt=dt+atmp(i,j)*xxx(j)
	    end do
	    disp1=disp1+abs(btmp(i))
	    disp2=disp2+abs(btmp(i)-dt)
	    !write(*,'(3f8.4)')bbb(i),bbb(i)-dt,dt
        end do
        disp1=disp1/(nk-3)
        disp2=disp2/(nk-3)
        red=100*(disp1-disp2)/disp1

        shift0=sqrt(xxx(1)*xxx(1)+xxx(2)*xxx(2))

        scale=1.
        if(shift0.gt.step_cur)scale=step_cur/shift0

        dxcur=-xxx(1)*scale
        dycur=-xxx(2)*scale
        dzcur=0

        xnew=xzt+dxcur
        ynew=yzt+dycur
        znew=zzt+dzcur
		
    end if
	!stop

    dscur=sqrt(dxcur*dxcur+dycur*dycur+dzcur*dzcur)

    xzt=xnew
    yzt=ynew
    zzt=znew

    if(key_ft1_xy2.eq.1) then
        call decsf(xzt,yzt,zzt,fi0,tet0,fff,ttt,hhh)
    else
        fff=xxt; ttt=yzt
    end if

    depth = sph_flat(key_flat1,xzt,yzt,zzt)

    limit=0
    zlim_up=h_lim(fff,ttt)
    !write(*,*)' depth=',depth,' zlim_up=',zlim_up
    if (depth.le.zlim_up) then
        depth=zlim_up
        zzt = flat_sph(key_flat1,xzt,yzt,depth)
        limit=1
    end if

    deplim = z_lim(fff,ttt)
    !write(*,*)' depth=',depth,' deplim=',deplim
    if (depth.ge.deplim) then
        depth=deplim
        zzt = flat_sph(key_flat1,xzt,yzt,depth)
        limit=1
    end if
    !write(*,*)' limit=',limit


    dstot=dstot+dscur

    depth = sph_flat(key_flat1,x_best,y_best,z_best)


    !if(iter.eq.6) stop
   ! write(*,*)' dscur=',dscur,' stepmin=',stepmin


    if(dscur.gt.stepmin) goto 331

    if(mod(nzt,ifreq).eq.0)write(*,*)x_best,y_best,z_best,depth
    if(key_true1.eq.1) write(11,*)xtrue,ytrue,ztrue
    write(11,*)x_best,y_best,depth,nkrat
    write(41,*)x_best,y_best,depth,nkrat
    !write(*,*)x_best,y_best,depth,nkrat


    if(key_ft1_xy2.eq.1) then
        call decsf(x_best,y_best,0.,fi0,tet0,fzt,tzt,hhh)
    else
        fzt=x_best; tzt=y_best
    end if

    write(14,*) fzt,tzt,z_best

    if(key_info.eq.1) then
        sec2=sec1-aver3
!        if(abs(aver).gt.10) then
!            write(*,*)' nztgood=',nztgood,' aver=',aver
!            write(*,*)' initial:',xold,yold,zold
!            write(*,*)' new:',xzt,yzt,zzt
!             do i=1,nkrat
!                write(*,*)ipskr(i),istkr(i),tobkr(i),trfkr(i)
!                dtot=dtot+abs(tobkr(i)-trfkr(i))
!                nray=nray+1
!                if(ipskr(i).eq.1) nrp=nrp+1
!                if(ipskr(i).eq.2) nrs=nrs+1
!            end do
!                   pause
!        end if
        write(15,'(i4,4(1x,i2),1x,2f8.2)')myr1,mmnt1,mdy1,mhr1,mmn1,sec2,sec1
        !call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
        write(16,*)fzt,tzt,zzt,nkrat
        write(16,'(i4,4(1x,i2),1x,2f8.2)')myr1,mmnt1,mdy1,mhr1,mmn1,sec2,sec1
        do i=1,nkrat
            write(16,*)stacod(istkr(i)),ipskr(i),tobkr(i),trfkr(i)
        end do

    end if


    disp1=0
    disp2=0
    do i=1,nkrat
        !write(*,*)istkr(i),ipskr(i),' dt=',tob_best(i)-trf_best(i),dtold(i)

        ist=istkr(i)
        ips=ipskr(i)
        !if(ist.ne.35.or.ips.ne.2) cycle
        xst=xstn(ist)
        yst=ystn(ist)
        zst=zstn(ist)
        dshor=sqrt((xst-xzt)*(xst-xzt)+(yst-yzt)*(yst-yzt))
        dsver=abs(zzt-zst)
        dist=sqrt(dshor*dshor+dsver*dsver)
        dt=tob_best(i)-trf_best(i)
        if(ips.eq.1) write(21,*)dshor,dt
        if(ips.eq.2) write(22,*)dshor,dt



        disp1=disp1+abs(dtold(i))
        disp2=disp2+abs(tob_best(i)-trf_best(i))
        dis_tot1=dis_tot1+abs(dtold(i))
        dis_tot2=dis_tot2+abs(tob_best(i)-trf_best(i))
        ntot=ntot+1
        write(11,*)ipskr(i),istkr(i),tob_best(i),trf_best(i)
        write(41,*)ipskr(i),istkr(i),tob_best(i)-trf_best(i)
        write(12)nod_best(i)
        if(ipskr(i).eq.1) nray_p=nray_p+1
        if(ipskr(i).eq.2) nray_s=nray_s+1
        !write(*,*)' nod_best(i)=',nod_best(i)
        do inod=1,nod_best(i)
            write(12)(rays_best(i3,inod,i),i3=1,3)
        end do
    end do
    disp1=disp1/nkrat
    disp2=disp2/nkrat
    s_best=sqrt((xini-x_best)**2+(yini-y_best)**2+(zini-z_best)**2)

    dt_aver1=dis_tot1/ntot
    dt_aver2=dis_tot2/ntot
    red=100*(dis_tot1-dis_tot2)/dis_tot1
    if(mod(nzt,ifreq).eq.0)write(*,*)' old resid=',dt_aver1,' new_resid=',dt_aver2,' red=',red
    if(mod(nzt,ifreq).eq.0)write(*,*)nzt,ntot,' ds=',s_best,' G=',goal_best
    if(mod(nzt,ifreq).eq.0)write(*,*)' ********************************************'

    deallocate(rays_all,rays_best,nod_all,nod_best)

        
    if(key_true1.eq.1) then
        err = sqrt ( (xtrue-x_best)**2 + (ytrue-y_best)**2 + (ztrue-depth)**2 )
        err_loc = err_loc + err
        err0 = sqrt ( (xtrue-xini)**2 + (ytrue-yini)**2 + (ztrue-depth0)**2 )
        err_loc0 = err_loc0 + err0
    end if


!if(nzt.eq.3) stop

goto 21
22 close(1)

close(11)
close(41)
close(12)
close(14)
close(21)
close(22)
if(key_info.eq.1)close(3)
if(key_info.eq.1)close(15)
if(key_info.eq.1)close(16)

open(11,file='../../../DATA/'//ar//'/'//md//'/data/numray'//itt//'.dat')
write(11,*) nray_p,nray_s
write(*,*) ' number of rays: P=',nray_p,' S=',nray_s
close(11)


if(key_true1.eq.1) then
    err_loc = err_loc / nzt
    err_loc0 = err_loc0 / nzt
    write(*,*)' err_loc0=',err_loc0,' err_loc=',err_loc
end if


write(*,*)' nzt=',nzt,' nray=',nray



key_preview=0
open(1,file='../../../preview_key.txt')
read(1,*,end=771)key_preview
771 close(1)



!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

if(key_preview.eq.0) stop

1313 continue
!call visual_srces_stat(ar,md,itt)


i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\LOC')


open(1,file='../../../DATA/'//ar//'/config.txt')
read(1,*)
read(1,*) npix_h_x0,npix_h_y0
read(1,*)tick_h_x,tick_h_y
read(1,*)
read(1,*) npix_v_x0,npix_v_z0
read(1,*)tick_v_x,tick_v_z
close(1)

open(1,file='../../../DATA/'//ar//'/sethor.dat')
read(1,*)
read(1,*)
read(1,*)xmap1,xmap2,dx,ymap1,ymap2,dy
close(1)

dxmark=0; dymark=0
open(1,file='../../../DATA/'//ar//'/setver.dat')
read(1,*) nvert
do ivert=1,nvert
    read(1,*)xvert1(ivert),yvert1(ivert),xvert2(ivert),yvert2(ivert)
end do
read(1,*)dist_max
read(1,*)dxsec0
read(1,*)zmap1,zmap2
read(1,*)dsmark
read(1,*)
read(1,*,end=991,err=991)dxmark,dymark
991 close(1)


i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')


!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
! 1: Vertical sections

do ivert=1,nvert


    write(ver,'(i2)')ivert

    if(key_ft1_xy2.eq.1) then
        fia=xvert1(ivert)
        fib=xvert2(ivert)
        teta=yvert1(ivert)
        tetb=yvert2(ivert)
            !write(*,*)' fia=',fia,' teta=',teta
            !write(*,*)' fib=',fib,' tetb=',tetb

        call SFDEC(fia,teta,0.,xa,ya,Z,fi0,tet0)
        call SFDEC(fib,tetb,0.,xb,yb,Z,fi0,tet0)

    else

        xa=xvert1(ivert)
        ya=yvert1(ivert)
        xb=xvert2(ivert)
        yb=yvert2(ivert)

    end if




    dist=sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya))
    sinpov=(yb-ya)/dist
    cospov=(xb-xa)/dist

    !write(*,*)' xa=',xa,' ya=',ya
    !write(*,*)' xb=',xb,' yb=',yb

    !write(*,*)' dist=',dist
    nxsec=dist/dxsec0+1
    dxsec=dist/(nxsec-1)

    !Draw topography on the section
    open(11,file='../../../TMP_files/loc/topo_'//ver//'.bln')
    write(11,*)nxsec
    do ix=1,nxsec
        sss=(ix-1)*dxsec
        xcur=xa+((xb-xa)/dist)*sss
        ycur=ya+((yb-ya)/dist)*sss
        if(key_ft1_xy2.eq.1) then
            call decsf(xcur,ycur,0.,fi0,tet0,fff,ttt,h)
        else
            fff=xcur
            ttt=ycur
        end if
        depth=h_lim(fff,ttt)
        ztopo=flat_sph(key_flat1,xcur,ycur,depth)
        write(11,*)sss,-ztopo
    end do
    close(11)



    open(11,file='../../../TMP_files/loc/mark_'//ver//'.dat')
    imark=0
    do sss=0.,dist,dsmark
        x=xa+cospov*sss
        y=ya+sinpov*sss
 
        if(key_ft1_xy2.eq.1) then
            call decsf(x,y,0,fi0,tet0,fi,tet,hhh)
            x=fi; y=tet
        end if
 
        write(11,*)x,y,sss
        imark=imark+1
        xmark(ivert,imark)=x
        ymark(ivert,imark)=y
        smark(ivert,imark)=sss
    end do
    imark=imark+1
        if(key_ft1_xy2.eq.1) then
            xmark(ivert,imark)=fib
            ymark(ivert,imark)=tetb
        else
            xmark(ivert,imark)=xb
            ymark(ivert,imark)=yb
       end if
    smark(ivert,imark)=dist
    close(11)
    nmark(ivert)=imark

	
    ! Draw the position of the section on the surface (line)
    open(11,file='../../../TMP_files/loc/mark_'//ver//'.bln')
    !write(*,*)' nmark(ivert)=',nmark(ivert)
    write(11,*) imark
    do i=1,imark

        write(11,*)xmark(ivert,i),ymark(ivert,i)	!,smark(i)
        !write(*,*)xmark(ivert,i),ymark(ivert,i)	!,smark(i)
    end do
    close(11)



    if(npix_v_x0.ne.0) then
        npix_x=npix_v_x0
    else
        npix_x=int(npix_v_z0*(dist/(zmap2-zmap1)))
        npix_z=npix_v_z0
    end if

    write(*,*)' Vertical:',ver,' dist=',dist,' npix_x=',npix_x,' npix_z=',npix_z

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
    open(11,file='../../../TMP_files/loc/stat_ver'//ver//'.dat')
    201	continue
        read(1,*,end=202)xst,yst,depth
        zst = flat_sph(key_flat1,xst,yst,depth)

        xx1=(xst-xa)*cospov+(yst-ya)*sinpov
        yy1=-(xst-xa)*sinpov+(yst-ya)*cospov

        if(abs(yy1).lt.dist_max) write(11,*)xx1,-zst,yy1
        goto 201
    202 close(1)
    close(11)

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays'//itt//'.dat')
    open(12,file='../../../TMP_files/loc/srces_ver'//itt//ver//'.dat')
    if(key_true1.eq.1) open(13,file='../../../TMP_files/loc/shift_ver'//itt//ver//'.bln')
    203	continue
        if(key_true1.eq.1)read(1,*,end=204)xzt0,yzt0,zzt0
        read(1,*,end=204)xzt,yzt,depth,nkrat
        zzt = flat_sph(key_flat1,xzt,yzt,depth)



        if(key_true1.eq.1) xx0=(xzt0-xa)*cospov+(yzt0-ya)*sinpov
        if(key_true1.eq.1) yy0=-(xzt0-xa)*sinpov+(yzt0-ya)*cospov
        xx1=(xzt-xa)*cospov+(yzt-ya)*sinpov
        yy1=-(xzt-xa)*sinpov+(yzt-ya)*cospov

 
        do i=1,nkrat
            read(1,*)ips,ist,tobs,tref
        end do
        if(abs(yy1).lt.dist_max) write(12,*)xx1,-zzt

        if(key_true1.eq.1) then
            if(abs(yy1).lt.dist_max) then
                write(13,*)2
                write(13,*)xx0,-zzt0
                write(13,*)xx1,-zzt
            end if
        end if

        goto 203
    204 close(1)
    close(12)
    if(key_true1.eq.1) close(13)



    open(14,file='config.txt')
    write(14,*)npix_x,npix_z
    write(14,*)'_______ Size of the picture in pixels (nx,ny)'
    write(14,*)0,dist
    write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
    write(14,*)-zmap2,-zmap1
    write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
    write(14,*)tick_v_x,tick_v_z
    write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
    write(14,205)ar,md,itt,ver
    205 format('..\..\..\PICS\',a8,'\',a8,'\LOC\loc_ver',a1,a2,'.png')
    write(14,*)'_______ Path of the output picture'
    if(key_true1.eq.1) then
        write(14,206)itt,ver,err_loc
        206 format(' Source locations in iteration: ',a1,' , vertical section',a2,' (err=',f7.2,' m)')
    else
        write(14,276)itt,ver
        276 format(' Source locations in iteration: ',a1,' , vertical section',a2)
    end if
    write(14,*)'_______ Title of the plot on the upper axe'
    write(14,*)	1
    write(14,*)'_______ Number of layers'


    write(14,59)
    write(14,*)	3
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,207)itt,ver
    207 format('..\..\..\TMP_files\loc\srces_ver',a1,a2,'.dat')
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
    write(14,208)ver
    208 format('..\..\..\TMP_files\loc\stat_ver',a2,'.dat')
    write(14,*)'_______ Location of the DAT file'
    write(14,*)	2
    write(14,*)'_______ Symbol (1: circle, 2: square)'
    write(14,*)	5
    write(14,*)'_______ Size of dots in pixels'
    write(14,*)	0,0,255
    write(14,*)'_______ RGB color'

    write(14,59)
    write(14,*)	2
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,214)ver
    214 format('..\..\..\TMP_files\loc\topo_',a2,'.bln')
    write(14,*)'_______ Location of the BLN file'
    write(14,*)	3
    write(14,*)'_______ Thickness of line in pixels'
    write(14,*)	0,0,0
    write(14,*)'_______ RGB color'


    if(key_true1.eq.1) then

        write(14,59)
        write(14,*)	2
        write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
        write(14,210)ver
        210 format('..\..\..\TMP_files\loc\shift_ver',a2,'.bln')
        write(14,*)'_______ Location of the BLN file'
        write(14,*)	1
        write(14,*)'_______ Thickness of line in pixels'
        write(14,*)	0,0,0
        write(14,*)'_______ RGB color'
    end if

    close(14)


    i=system('layers.exe')


end do !vertical sections


!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************

! 1: Horizontal section

if(npix_h_x0.ne.0) then
    npix_x=npix_h_x0
else
    npix_x=int(npix_h_y0*((xmap2-xmap1)/(ymap2-ymap1)))
    npix_y=npix_h_y0
end if

if(npix_h_y0.ne.0) then
    npix_y=npix_h_y0
else
    npix_y=int(npix_h_x0*((ymap2-ymap1)/(xmap2-xmap1)))
    npix_x=npix_h_x0
end if

write(*,*)' Horizontal: npix_x=',npix_x,' npix_y=',npix_y

if(key_ft1_xy2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/inidata/stat_ft.dat')
else
    open(1,file='../../../DATA/'//ar//'/inidata/stat_xy.dat')
end if
open(11,file='../../../TMP_files/loc/stat_hor.dat')
nst=0
431	continue
    read(1,*,end=432)xst,yst,zst
    write(11,*)xst,yst
    nst=nst+1
    goto 431
432 close(1)
close(11)
write(*,*)' nst=',nst

open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays'//itt//'.dat')
open(12,file='../../../TMP_files/loc/srces_hor'//itt//'.dat')
if(key_true1.eq.1) open(13,file='../../../TMP_files/loc/shift_hor'//itt//'.bln')
433	continue
    if(key_true1.eq.1) read(1,*,end=434)xzt0,yzt0,zzt0
    read(1,*,end=434)xzt,yzt,zzt,nkrat
    do i=1,nkrat
        read(1,*)ips,ist,tobs,tref
    end do
    if(key_ft1_xy2.eq.1) then
      call decsf(xzt,yzt,zzt,fi0,tet0,fzt,tzt,hhh)
      call decsf(xzt0,yzt0,zzt,fi0,tet0,fzt0,tzt0,hhh)
      write(12,*)fzt,tzt
        if(key_true1.eq.1) then
            write(13,*)2
            write(13,*)fzt0,tzt0
            write(13,*)fzt,tzt
        end if
    else
        write(12,*)xzt,yzt
        if(key_true1.eq.1) then
            write(13,*)2
            write(13,*)xzt0,yzt0
            write(13,*)xzt,yzt
        end if
    end if
    goto 433
434 close(1)
close(12)
if(key_true1.eq.1) close(13)


open(14,file='config.txt')
write(14,*)npix_x,npix_y
write(14,*)'_______ Size of the picture in pixels (nx,ny)'
write(14,*)xmap1+fi0,xmap2+fi0
write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
write(14,*)ymap1+tet0,ymap2+tet0
write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
write(14,*)tick_h_x,tick_h_y
write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
write(14,151)ar,md,itt
151 format('..\..\..\PICS\',a8,'\',a8,'\LOC\loc_hor',a1,'.png')
write(14,*)'_______ Path of the output picture'
if(key_true1.eq.1) then
    write(14,611)itt,err_loc
    611 format(' Source locations in iteration: ',a1,' , map view (err=',f7.2,' m)')
else
    write(14,613)itt
    613 format(' Source locations in iteration: ',a1,' , map view ')
end if
write(14,*)'_______ Title of the plot on the upper axe'
write(14,*)	1
write(14,*)'_______ Number of layers'

59 format('********************************************')

write(14,59)
write(14,*)	3
write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
write(14,51)itt
51 format('..\..\..\TMP_files\loc\srces_hor',a1,'.dat')
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
write(14,52)
52 format('..\..\..\TMP_files\loc\stat_hor.dat')
write(14,*)'_______ Location of the DAT file'
write(14,*)	2
write(14,*)'_______ Symbol (1: circle, 2: square)'
write(14,*)	5
write(14,*)'_______ Size of dots in pixels'
write(14,*)	0,0,255
write(14,*)'_______ RGB color'

if(key_true1.eq.1) then
    write(14,59)
    write(14,*)	2
    write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
    write(14,621)itt
    621 format('..\..\..\TMP_files\loc\shift_hor',a1,'.bln')
    write(14,*)'_______ Location of the BLN file'
    write(14,*)	1
    write(14,*)'_______ Thickness of line in pixels'
    write(14,*)	0,0,0
    write(14,*)'_______ RGB color'
end if

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


do iver=1,nvert
	write(ver,'(i2)')iver

	write(14,59)
	write(14,*)	2
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,154)ver
154 format('..\..\..\TMP_files\LOC\mark_',a2,'.bln')
	write(14,*)'_______ Location of the BLN file'
	write(14,*)	2
	write(14,*)'_______ Thickness of line in pixels'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	write(14,59)
	write(14,*)	3
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,156)ver
156 format('..\..\..\TMP_files\LOC\mark_',a2,'.dat')
	write(14,*)'_______ Location of the DAT file'
	write(14,*)	1
	write(14,*)'_______ Symbol (1: circle, 2: square)'
	write(14,*)	6
	write(14,*)'_______ Size of dots in pixels'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	xa=xmark(iver,1)+dxmark
	ya=ymark(iver,1)+dymark
 
	write(14,59)
	write(14,*)	5
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,756)iver
756 format(i2,' - A')
	write(14,*)'_______ Text to output'
	write(14,*)	xa,ya
	write(14,*)'_______ Coordinates'
	write(14,*)	1
	write(14,*)'_______ Type (1: Arial, 2: Times New Romans)'
	write(14,*)	12
	write(14,*)'_______ Size'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	xb=xmark(iver,nmark(iver))+dxmark
	yb=ymark(iver,nmark(iver))+dymark

	write(14,59)
	write(14,*)	5
	write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
	write(14,757)iver
757 format(i2,' - B')
	write(14,*)'_______ Text to output'
	write(14,*)	xb,yb
	write(14,*)'_______ Coordinates'
	write(14,*)	1
	write(14,*)'_______ Type (1: Arial, 2: Times New Romans)'
	write(14,*)	12
	write(14,*)'_______ Size'
	write(14,*)	0,0,0
	write(14,*)'_______ RGB color'

	do imark=2,nmark(iver)-2

		xx=xmark(iver,imark)+dxmark
		yy=ymark(iver,imark)+dymark
		isss=smark(iver,imark)

		write(14,59)
		write(14,*)	5
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,'(i4)')isss
		write(14,*)'_______ Text to output'
		write(14,*)	xx,yy
		write(14,*)'_______ Coordinates'
		write(14,*)	1
		write(14,*)'_______ Type (1: Arial, 2: Times New Romans)'
		write(14,*)	10
		write(14,*)'_______ Size'
		write(14,*)	0,0,0
		write(14,*)'_______ RGB color'

	end do

end do


close(14)


i=system('layers.exe')





stop
end