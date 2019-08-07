subroutine	goal_new(xzt,yzt,zzt, disp,aver,nk,ank)

real tall(20),hall(20),aall(20)
real trfkr1(500),alkr1(500)
real dtk(500),dtk1(500),dtall(20,500),trfal(20,500),alall(20,500),dtmp(500)
integer kmin(500),nall(500),ngmin(500)

common/krat/nkrat,istkr(500),tobkr(500),ipskr(500),qualkr(500),trfkr(500),ngood(500),alkr(500),diskr(500)
common/iprint/iprint
common/center/f0,t0
common/stations/ xstat(9000),ystat(9000),zstat(9000)
common/stat_level/stat_level

common/loc_param/wgs,res_loc1,res_loc2,dist_limit,n_pwr_dist,ncyc_av,w_P_S_diff


REAL PI/3.1415926/
PER=PI/180.
rz=6371.

!write(*,*)' resmax=',res_max_loc


!dhzt=dh_crust_sm(xzt,yzt)
!write(*,*)' dhzt=',dhzt,' zzt=',zzt
!pause
!if(dhzt.ne.0.)pause

!write(*,*)' zzt=',zzt

vzt1=vrefmod(zzt,1)
vzt2=vrefmod(zzt,2)
vst1=vrefmod(stat_level,1)
vst2=vrefmod(stat_level,2)
!write(*,*)' vst1=',vst1,' vst2=',vst2
rzt=rz-zzt
rst=rz

do ikr=1,nkrat
    ips=ipskr(ikr)
    ist=istkr(ikr)
    !write(*,*)' ikr=',ikr,' ips=',ips,' ist=',ist
    xst=xstat(ist)
    yst=ystat(ist)
    zst=zstat(ist)
    vzt=vzt1
    vst=vst1
    if(ips.eq.2) then
	    vzt=vzt2
	    vst=vst2
    end if

    dist=sqrt((xst-xzt)*(xst-xzt)+(yst-yzt)*(yst-yzt))
    diskr(ikr)=sqrt(dist*dist+zzt*zzt)

    if(diskr(ikr).lt.0.1) then
	    tttt=0
	    aaaa=90
	    goto 441
    end if
    !write(*,*)

    !write(*,*)' dis=',dist,' zzt=',zzt

    call refmod_all(dist,zzt,ips, nal,tall,hall,aall)

    !write(*,*)' tob=',tobkr(ikr),' trf=',(tall(i),i=1,nal)
    !pause


    ! Select the first arrival
    tmin=999999
    do ii=1,nal
	    if(tall(ii).gt.tmin) cycle
	    tmin=tall(ii)
	    imin=ii
    end do
    hhhh=hall(imin)
    aaaa=aall(imin)
    tttt=tmin

    441	continue

    px=sin(aaaa*per)/vzt
    alkr(ikr)=aaaa
    !write(*,*)' px=',px,' aaaa=',aaaa

! Correction for the station elevation	
    slowst=1/vst
    sqr=slowst*slowst-px*px

    dzstat = zst - stat_level

    if(sqr.le.0) then
	    dtstat=0
    else
	    dtstat=-dzstat*sqrt(sqr)
    end if
    !write(*,*)' zst=',zst,vst,' px=',px,' dtstat=',dtstat
    !pause
    tttt = tttt + dtstat

    trfkr(ikr)=tttt
    alkr(ikr)=aaaa
    !write(*,*)' dis=',diskr(ikr),ipskr(ikr),' tob=',tobkr(ikr),' trf=',trfkr(ikr)
    dtk(ikr)=(tobkr(ikr)-tttt)

end do

do ikr=1,nkrat
	if(ipskr(ikr).eq.1) cycle
! 3 means that a differential time Tp-Ts is used
	do ik2=1,nkrat
		if(ipskr(ik2).eq.2) cycle
		if(istkr(ik2).ne.istkr(ikr)) cycle
		ipskr(ikr) = 3
		ikr2=ik2
		dtk(ikr)=dtk(ikr)-dtk(ikr2)
		goto 771
	end do
771	continue
end do

!do ikr=1,nkrat
	!write(*,*)' ips=',istkr(ikr),ipskr(ikr),' dis=',diskr(ikr),' dt=',dtk(ikr)
!end do
call dispers(dtk, disp,aver,nk,ank)
!write(*,*)' ank=',ank
!pause


do i=1,nkrat
	if(ipskr(i).eq.3) ipskr(i)=2
end do

return
end