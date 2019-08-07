subroutine dispers(dtk,	disp,aver,nk,goal)

real dtk(500)
integer ngood2(500)
common/loc_param/wgs,res_loc1,res_loc2,dist_limit,n_pwr_dist,ncyc_av,w_P_S_diff
common/krat/nkrat,istkr(500),tobkr(500),ipskr(500),qualkr(500),trfkr(500),ngood(500),alkr(500),diskr(500)
common/iprint/iprint


!n_pwr_dist=1
!ncyc_av=10
!w_P_S_diff=2

ngood=1
! Average residual :
313 continue
aver=0
nav=0
aver=0
nav=0
do i=1,nkrat
    if(ipskr(i).eq.3) cycle
    if(ngood(i).eq.0) cycle
    aver=aver+dtk(i)
    nav=nav+1
end do
aver=aver/nav
!write(*,*)' aver=',aver,' nav=',nav

!if(nav.gt.2) then
    dmax=-9999
    do i=1,nkrat
        if(ipskr(i).eq.3) cycle
        if(ngood(i).eq.0) cycle
        dtt=abs(dtk(i)-aver)
        !write(*,*)' i=',i,' dtt=',dtk(i)-aver
        if(dtt.gt.dmax) then
            dmax=dtt
            imax=i
        end if
    end do
    !write(*,*)' dmax=',dmax,' imax=',imax,' res_loc2=',res_loc2
    !pause
    if(dmax.gt.res_loc2) then
        ngood(imax)=0
        goto 313
    end if
!else
!    return
!end if

goal=0
www=0
do i=1,nkrat

    dt=abs(dtk(i)-aver)
    if(ipskr(i).eq.3)dt=abs(dtk(i)) 

    if(dt.gt.res_loc2) then
	    aaa=0
    else if(dt.lt.res_loc1) then
	    aaa=1
    else
	    aaa=(dt-res_loc2)/(res_loc1-res_loc2)
    end if

    ccc=1
    if(ipskr(i).eq.3) ccc = w_P_S_diff


    goal = goal + aaa*ccc
    www = www + ccc
    !write(*,'(2i4,5f9.3)')istkr(i),ipskr(i),dtk(i)-aver,aaa,ccc,goal
end do
goal = goal / www
!write(*,*)' goal=',goal
!pause




! Compute the dispersion :
disp=0
nk=0
do i=1,nkrat
	if(ngood(i).eq.0) cycle
	nk=nk+1
	dt=dtk(i)
	if(ipskr(i).ne.3) dt=(dtk(i)-aver)
	disp=disp+abs(dt)
	!write(*,'(2i4,5f9.3)')istkr(i),ipskr(i),dt
end do
disp=disp/nk
!write(*,*)' disp=',disp

return
end