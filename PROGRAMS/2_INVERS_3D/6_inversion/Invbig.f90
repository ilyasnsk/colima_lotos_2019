character*2 itinv
character*8 ar,md,line
character*1 it,itold,ppss,rm,gr

real matruzel(3000)
integer muzel(3000),ist_p(5000),ist_s(5000)

real dv_zn_p(3),dv_zn_s(3)

real ylev_p(200),ylev_s(200)
integer nt_p(200),npz_p(10),nt_s(200),npz_s(10),nrps(2)

allocatable  xtop_p(:,:), ztop_p(:,:)
allocatable  xtop_s(:,:), ztop_s(:,:)
allocatable kpop_p(:,:),kobr_p(:,:),kpz_p(:)
allocatable kpop_s(:,:),kobr_s(:,:),kpz_s(:)

allocatable dvp_it(:),dvs_it(:),dt_it(:)
allocatable x_ztr(:),y_ztr(:),z_ztr(:),t_ztr(:),p_sta(:),s_sta(:)

allocatable x(:),x0(:),u(:),v(:),w(:),aaa(:),dt(:),xmod(:)
allocatable ncolrow(:),ncol(:)

common/velmod/nnps,ypros,xpros,zpros,cfps
common/inipar/nzone,nsurf,smth,xlim1,xlim2,ylim1,ylim2,zlim1,zlim2
common/itstep/itstep


pi=3.1415926
iter_max=1
it_curr=1
w_act=1


open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)iter		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(gr,'(i1)')igr
write(it,'(i1)')iter
write(itold,'(i1)')iter-1

write(*,*)' execution of invers'
write(*,*)' ar=',ar,' md=',md,' it=',it,' gr=',gr

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
read(1,*,end=441,err=441)key_ft1_xy2
read(1,*,end=441,err=441)key_true1
read(1,*,end=441,err=441)key_flat1      ! 1: calculations in flat model, 2: spherical velocity
441 close(1)


!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=573)line
	if(line.eq.'INVERSIO') goto 574
end do
573 continue
write(*,*)' cannot find INVERSION PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
574 continue
read(1,*)iter_lsqr
read(1,*)wg_vel_p,wg_vel_s
read(1,*)sm_hor_p,sm_hor_s
read(1,*)sm_ver_p,sm_ver_s
read(1,*)rg_amp_p,rg_amp_s
read(1,*)
read(1,*)wg_st_p,wg_st_s
read(1,*)wzt_hor
read(1,*)wzt_ver
read(1,*)wzt_time
close(1)
!******************************************************************

call read_vref(ar,md)
call read_3D_mod_v(ar,md,iter-1)
call read_ini_model(ar,md)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray1.dat')
read(1,*) nrps(1),nrps(2)
close(1)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numbers'//it//gr//'.dat')
read(1,*)nray,nobr_p,nobr_s,nztr,nonzer
read(1,*)nrp,nonz_p
read(1,*)nrs,nonz_s
close(1)

write(*,*)' Number of events =',nztr

allocate (x_ztr(nztr),y_ztr(nztr),z_ztr(nztr),t_ztr(nztr))
x_ztr=0
y_ztr=0
z_ztr=0
t_ztr=0


! Read the coordinates of the stations
open(2,file='../../../DATA/'//ar//'/'//md//'/data/stat_xy.dat')
i=0
3	i=i+1
	read(2,*,end=4)xstn,ystn,zstn
	goto 3
4	close(2)
nstat=i-1
write(*,*)' Number of stations=',nstat

allocate (p_sta(nstat),s_sta(nstat))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nrps(1).ne.0) then
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr1'//gr//'.dat')
	read(1,*) nvel_p
	allocate(dvp_it(nvel_p),kpz_p(nvel_p),kobr_p(nvel_p,2))
	read(1,*)((kobr_p(i,j),i=1,nvel_p), j=1,2)
	close(1)
end if

if(nrps(2).ne.0) then
	open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr2'//gr//'.dat')
	read(1,*) nvel_s
	allocate(dvs_it(nvel_s),kpz_s(nvel_s),kobr_s(nvel_s,2))
	read(1,*)((kobr_s(i,j),i=1,nvel_s), j=1,2)
	close(1)
end if
write(*,*)' nvel_p=',nvel_p,' nvel_s=',nvel_s

dvp_it=0
dvs_it=0

open(1,file='../../../DATA/'//ar//'/'//md//'/data/kall_xy1'//gr//'.dat')
read(1,*)nxpl,nypl
close(1)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nrps(1).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr1'//gr//'.dat')
    nmax=0
    do n=1,nypl
        read(1,*)nt,ylev_p(iy)
        if(nt.gt.nmax) nmax=nt
        do i=1,nt
            read(1,*)x,z,l
        end do
    end do
    close(1)
    !write(*,*)' nmax=',nmax

    allocate(xtop_p(nmax,nypl),ztop_p(nmax,nypl),kpop_p(nmax,nypl))

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr1'//gr//'.dat')
    npz_p=0
    do n=1,nypl
        read(1,*)nt
        !write(*,*)n,' ntop=',nt
        nt_p(n)=nt
        do i=1,nt_p(n)
            read(1,*)xtop_p(i,n),ztop_p(i,n),kpop_p(i,n)
        end do
    end do
    close(1)
end if

if(nrps(2).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr2'//gr//'.dat')
    nmax=0
    do n=1,nypl
        read(1,*)nt
        if(nt.gt.nmax) nmax=nt
        do i=1,nt
            read(1,*)x,z,l
        end do
    end do
    close(1)
    !write(*,*)' nmax=',nmax

    allocate(xtop_s(nmax,nypl),ztop_s(nmax,nypl),kpop_s(nmax,nypl))

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr2'//gr//'.dat')
    do n=1,nypl
        read(1,*)nt
        nt_s(n)=nt
        do i=1,nt_s(n)
            read(1,*)xtop_s(i,n),ztop_s(i,n),kpop_s(i,n)
        end do
    end do
    close(1)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nonz=0
ir=0
avdt=0
totmat=0
nst_p=0
nst_p=0

open(3,file='../../../TMP_files/tmp/matr'//it//gr//'.dat',form='binary')

95  read(3,end=96) nuz,resid,ist,ips,izt
    !write(*,*) nuz,resid,ist,ips,izt
    read(3)dtdx,dtdy,dtdz
    if(nuz.ne.0)read(3) (muzel(ii),matruzel(ii),ii=1,nuz)

    ir=ir+1
    nonz = nonz + nuz 
    if(ist.eq.0) goto 95

    nonz = nonz + 5

    if(ips.eq.1)then
        if(nst_p.ne.0) then
            do i=1,nst_p
	            if(ist_p(i).eq.ist) goto 298
            end do
        end if
        nst_p=nst_p+1
        ist_p(nst_p)=ist
    else
        if(nst_s.ne.0) then
            do i=1,nst_s
	            if(ist_s(i).eq.ist) goto 298
            end do
        end if
        nst_s = nst_s + 1
        ist_s(nst_s) = ist
    end if

    298	continue
    !if(ir.eq.35375)write(*,*)' ir=',ir,' nonz=',nonz

    goto 95
96  close(3)

write(*,*)' Main matrix: nonz=',nonz,' ir=',ir

allocate (dt_it(ir))
write(*,*)' Number of rays =',nray,ir,' nzt=',nztr
write(*,*)' nst_p=',nst_p,' nst_s=',nst_s


dt_it=0
p_sta=0
s_sta=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nc = nvel_p + nvel_s + nst_p + nst_s + nztr*4

nonzer=nonz
!write(*,*)nvel_p,nvel_s,nst_p,nst_s,nztr
!write(*,*)' number of columns: =',nc


np1=nvel_p
np2=np1+nvel_s
np3=np2+nst_p
np4=np3+nst_s
!write(*,*)' np1=',np1,' np2=',np2,' np3=',np3

nrows=nray

if(nrps(1).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_hor1'//gr//'.dat',form='binary')
    read(1)link_hor_p
    write(*,*)' link_hor_p=',link_hor_p
    close(1)
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_ver1'//gr//'.dat',form='binary')
    read(1)link_ver_p
    write(*,*)' link_ver_p=',link_ver_p
    close(1)
end if

if(nrps(2).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_hor2'//gr//'.dat',form='binary')
    read(1)link_hor_s
    write(*,*)' link_hor_s=',link_hor_s
    close(1)
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_ver2'//gr//'.dat',form='binary')
    read(1)link_ver_s
    write(*,*)' link_ver_s=',link_ver_s
    close(1)
end if

if(sm_hor_p.gt.0.0001.and.nrps(1).ne.0) then
	nonzer = nonzer + link_hor_p*2
	nrows = nrows + link_hor_p
end if
if(sm_ver_p.gt.0.0001.and.nrps(1).ne.0) then
	nonzer = nonzer + link_ver_p*2
	nrows = nrows + link_ver_p
end if
if(sm_hor_s.gt.0.0001.and.nrps(2).ne.0) then
	nonzer = nonzer + link_hor_s*2
	nrows = nrows + link_hor_s
end if
if(sm_ver_s.gt.0.0001.and.nrps(2).ne.0) then
	nonzer = nonzer + link_ver_s*2
	nrows = nrows + link_ver_s
end if
write(*,*)' After Smoothing : ir=',nrows,' nonz=',nonzer

if(rg_amp_p.gt.0.0001.and.nrps(1).ne.0) then
	nonzer=nonzer+nvel_p
	nrows=nrows+nvel_p
end if
if(rg_amp_s.gt.0.0001.and.nrps(2).ne.0) then
	nonzer=nonzer+nvel_s
	nrows=nrows+nvel_s
end if

write(*,*)' nvel_p=',nvel_p,' nvel_s=',nvel_s
write(*,*)' Preliminary values: N rows=',nrows,' N nonzer=',nonzer
write(*,*)

allocate(xmod(nc),x(nc),x0(nc),v(nc),w(nc))
allocate(ncol(nrows),dt(nrows),u(nrows))
allocate(aaa(nonzer),ncolrow(nonzer))

nonz=0  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ir0=0
ir=0
amatsum=0
open(3,file='../../../TMP_files/tmp/matr'//it//gr//'.dat',form='binary')
    97	read(3,end=98) nuz,resid,ist,ips,izt
    read(3)dtdx,dtdy,dtdz
    ir0=ir0+1
    if(nuz.ne.0)read(3) (muzel(ii),matruzel(ii),ii=1,nuz)
    ir=ir+1

!    if(ir.eq.24310) then
!        write(*,*)' nuz=',nuz,' resid=',resid
!        write(*,*)' ist=',ist,' ips=',ips,' izt=',izt
!        write(*,*)' dtdx=',dtdx,' dtdy=',dtdy,' dtdz=',dtdz
!        do ii=1,nuz
!            write(*,*)muzel(ii),matruzel(ii)
!        end do
!        stop
!    end if

    dt(ir) = resid  
    do ii=1,nuz
        mmm=muzel(ii)
        nnp=0
        if(ips.eq.2) nnp=np1
        add = wg_vel_p
        if(ips.eq.2) add = wg_vel_s
        nonz=nonz+1
        aaa(nonz) = matruzel(ii) * add 
        ncolrow(nonz)=nnp+mmm
    end do
    !write(*,*)' dt=',dt(ir)

    ncol(ir) = nuz 

    if(ist.eq.0) goto 97

    ncol(ir) = ncol(ir) + 4 + 1 

    if(ips.eq.1) then
        do i=1,nst_p
            if(ist_p(i).ne.ist) cycle
            nonz=nonz+1
            aaa(nonz)=wg_st_p
            ncolrow(nonz)=np2+i
            exit
        end do
    else
        do i=1,nst_s
            if(ist_s(i).ne.ist) cycle
            nonz=nonz+1
            if(nonz.ge.nonzer) pause'7' 
            aaa(nonz)=wg_st_s
            ncolrow(nonz)=np3+i
            exit
        end do
    end if

    nonz=nonz+1
    aaa(nonz)=dtdx*wzt_hor	
    ncolrow(nonz)=np4 + (izt-1)*4 + 1
    !write(*,*)' ztr:',nonz,ncolrow(nonz),aaa(nonz)
    nonz=nonz+1
    aaa(nonz)=dtdy*wzt_hor	
    ncolrow(nonz)=np4 + (izt-1)*4 + 2
    !write(*,*)' ztr:',nonz,ncolrow(nonz),aaa(nonz)
    nonz=nonz+1
    aaa(nonz)=dtdz*wzt_ver	
    ncolrow(nonz)=np4 + (izt-1)*4 + 3
    !write(*,*)' ztr:',nonz,ncolrow(nonz),aaa(nonz)
    nonz=nonz+1
    aaa(nonz)=wzt_time	
    ncolrow(nonz)=np4 + (izt-1)*4 + 4
    !write(*,*)' ztr:',nonz,ncolrow(nonz),aaa(nonz)
    !pause


    goto 97
98    continue
close(3)
write(*,*)' Main matrix: nonz=',nonz,' ir=',ir

if(sm_hor_p.gt.0.0001.and.nrps(1).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_hor1'//gr//'.dat',form='binary')
    read(1)link_hor_p
    do ilink=1,link_hor_p
        read(1)ipar1,ipar2
        ir=ir+1
        dt(ir)=0.
        ncol(ir)=2
        smth=sm_hor_p
        nonz=nonz+1
        aaa(nonz)=smth; ncolrow(nonz) = ipar1
        nonz=nonz+1
        aaa(nonz)=-smth; ncolrow(nonz) = ipar2
    end do
    close(1)
end if
if(sm_ver_p.gt.0.0001.and.nrps(1).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_ver1'//gr//'.dat',form='binary')
    read(1)link_ver_p
    do ilink=1,link_ver_p
        read(1)ipar1,ipar2
        ir=ir+1
        dt(ir)=0.
        ncol(ir)=2
        smth=sm_ver_p
        nonz=nonz+1
        aaa(nonz)=smth; ncolrow(nonz) = ipar1
        nonz=nonz+1
        aaa(nonz)=-smth; ncolrow(nonz) = ipar2
    end do
    close(1)
end if
if(sm_hor_s.gt.0.0001.and.nrps(2).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_hor2'//gr//'.dat',form='binary')
    read(1)link_hor_s
    do ilink=1,link_hor_s
        read(1)ipar1,ipar2
        ir=ir+1
        dt(ir)=0.
        ncol(ir)=2
        smth=sm_hor_s
        nonz=nonz+1
        aaa(nonz)=smth; ncolrow(nonz) = np1 + ipar1
        nonz=nonz+1
        aaa(nonz)=-smth; ncolrow(nonz) = np1 + ipar2
    end do
    close(1)
end if
if(sm_ver_s.gt.0.0001.and.nrps(2).ne.0) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/data/link_ver2'//gr//'.dat',form='binary')
    read(1)link_ver_s
    do ilink=1,link_ver_s
        read(1)ipar1,ipar2
        ir=ir+1
        dt(ir)=0.
        ncol(ir)=2
        smth=sm_ver_s
        nonz=nonz+1
        aaa(nonz)=smth; ncolrow(nonz) = np1 + ipar1
        nonz=nonz+1
        aaa(nonz)=-smth; ncolrow(nonz) = np1 + ipar2
    end do
    close(1)
end if

write(*,*)' After Smoothing : ir=',ir,' nonz=',nonz

if(rg_amp_p.gt.0.0001.and.nrps(1).ne.0) then
    do i=1,nvel_p 
        ir=ir+1
        dt(ir)=0.
        ncol(ir)=1
        nonz=nonz+1
        aaa(nonz)=rg_amp_p
        ncolrow(nonz)= i
    end do
end if

if(rg_amp_s.gt.0.0001.and.nrps(2).ne.0) then
    do i=1,nvel_s 
        ir=ir+1
        dt(ir)=0.
        ncol(ir)=1
        nonz=nonz+1
        aaa(nonz)=rg_amp_s
        ncolrow(nonz)= np1 + i
    end do
end if

write(*,*)' After Regularization: ir=',ir,' nonz=',nonz

nz=nonz
nr=ir

do i=1,nr
	u(i)=dt(i)
end do


kount=0
amatsum=0
do irw=1,nr
    do ii=1,ncol(irw)
        kount=kount+1
        nn=ncolrow(kount)
        !write(*,*)' nn=',nn,' aaa=',aaa(kount)
        if(abs(u(irw)).gt.1.e10) then
            write(*,*)' Warning!!!'
            write(*,*)' dt=',u(irw)
            pause
        end if
        amatsum=amatsum+aaa(kount)
        if(abs(aaa(kount)).gt.1.e10) then
            write(*,*)' Warning!!!'
            write(*,*)' irw=',irw,' kount=',kount,' nn=',nn
            write(*,*)' aaa=',aaa(kount)
            stop
        end if
        if(nn.le.0.or.nn.gt.nc) then
            write(*,*)' Warning!!!'
            write(*,*)' irw=',irw,' nn=',nn
            pause
        end if
    end do
    !if(irw.ge.31400.and.irw.le.31500) write(*,*)' irw=',irw,' dt=',u(irw),amatsum
end do
write(*,*)' N rows=',nr,' N columns=',nc,' N nonzer=',nz
write(*,*)

!*******************************************************
!*******************************************************
!*******************************************************
call pstomo(nr,nc,x,u,v,w,iter_lsqr,nz,aaa,ncolrow,ncol)
!*******************************************************
!*******************************************************
!*******************************************************
err=0
kount=0
avres=0.
avres0=0.
dtvel=0
dtsta=0
dtztr=0
write(*,*)' nray=',nray
do irw=1,nray
	dt1=0
	do ii=1,ncol(irw)
		kount=kount+1
		nn=ncolrow(kount)
		dt1=dt1+aaa(kount)*x(nn)
	end do
	if(irw.le.nray)dt_it(irw)=dt_it(irw)+dt1
	avres=avres+abs(dt(irw)-dt1)
	avres0=avres0+abs(dt(irw))
end do
avres=avres/nray
avres0=avres0/nray
reduct=((avres0-avres)/avres0)*100.

write(*,*)'_____________________________________________________________'
write(*,*)' avres0=',avres0,' avres=',avres,' red=',reduct
write(*,*)'_____________________________________________________________'


if(nrps(1).ne.0) then 
	avdvp=0
	open(11,file='../../../DATA/'//ar//'/'//md//'/data/vel_p_'//it//gr//'.dat')
	if(iter.gt.1) open(1,file='../../../DATA/'//ar//'/'//md//'/data/vel_p_'//itold//gr//'.dat')
	do i=1,nvel_p

		if(iter.gt.1) then
			read(1,*)dvp_old,vvp
		else
			dvp_old=0
			i1=kobr_p(i,1)
			i2=kobr_p(i,2)
			!write(*,*)' i1=',i1,' i2=',i2
			if(i1.eq.0) cycle
			xx=xtop_p(i1,i2)
			yy=ylev_p(i2)
			zz=ztop_p(i1,i2)

			vvp=velocity(xx,yy,zz,1)
			!write(*,*)' xx=',xx,' yy=',yy,' zz=',zz
			!write(*,*)' vvp=',vvp
		end if

		dvp = x(i) * wg_vel_p

		dvp_it(i) = dvp_it(i) + dvp

		dvp_new = dvp_old + dvp_it(i)

		write(11,*)dvp_new,vvp

		avdvp=avdvp+100*abs(dvp)/vvp
	end do
	close(11)
	if(iter.gt.1) close(1)
	avdvp=avdvp/nvel_p
end if

!write(*,*)' nrps=',nrps(1),nrps(2)

if(nrps(2).ne.0) then 
	avdvs=0
	open(11,file='../../../DATA/'//ar//'/'//md//'/data/vel_s_'//it//gr//'.dat')
	if(iter.gt.1) open(1,file='../../../DATA/'//ar//'/'//md//'/data/vel_s_'//itold//gr//'.dat')
	do i=1,nvel_s

		if(iter.gt.1) then
			read(1,*)dvs_old,vvs
		else
			dvs_old=0
			i1=kobr_s(i,1)
			i2=kobr_s(i,2)
			!write(*,*)' i1=',i1,' i2=',i2
			if(i1.eq.0) cycle
			xx=xtop_s(i1,i2)
			yy=ylev_s(i2)
			zz=ztop_s(i1,i2)

			!write(*,*)' xx=',xx,' yy=',yy,' zz=',zz
			vvs=velocity(xx,yy,zz,2)
		end if

		dvs = x(np1+i) * wg_vel_s

		dvs_it(i) = dvs_it(i) + dvs

		dvs_new = dvs_old + dvs_it(i)

		write(11,*)dvs_new,vvs

		avdvs=avdvs+100*abs(dvs)/vvs
	end do
	close(11)
	if(iter.gt.1) close(1)

	avdvs=avdvs/nvel_s
	avdv=0
	do i=1,nst_s
		ii=np3+i
		avdv=avdv+abs(x(ii)*wg_st_s)
	end do
	avstat_s=avdv/nst_s
	open(12,file='../../../DATA/'//ar//'/'//md//'/data/stcor_s_'//it//gr//'.dat')
	do ist=1,nstat
		s_sta(ist)=0
		do i=1,nst_s
			if(ist_s(i).ne.ist) cycle
			ii=np3+i
			s_sta(ist)=x(ii)*wg_st_s
			exit
		end do
		write(12,*)s_sta(ist)
	end do
	close(12)
end if

write(*,*)' total avdv_p=',avdvp,' avdv_s=',avdvs


avdv=0
do i=1,nst_p
	ii=np2+i
	avdv=avdv+abs(x(ii)*wg_st_p)
end do
avstat_p=avdv/nst_p


open(12,file='../../../DATA/'//ar//'/'//md//'/data/stcor_p_'//it//gr//'.dat')
do ist=1,nstat
	p_sta(ist)=0
	do i=1,nst_p
		if(ist_p(i).ne.ist) cycle
		ii=np2+i
		p_sta(ist)=x(ii)*wg_st_p
	end do
	write(12,*)p_sta(ist)
end do
close(12)

write(*,*)' avstat_p=',avstat_p,' avstat_s=',avstat_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


open(12,file='../../../DATA/'//ar//'/'//md//'/data/ztcor_'//it//gr//'.dat')
av1=0
av2=0
av3=0
do i=1,nztr
	ii=np4+(i-1)*4
	x_ztr(i)=x_ztr(i)+x(ii+1)*wzt_hor
	y_ztr(i)=y_ztr(i)+x(ii+2)*wzt_hor
	z_ztr(i)=z_ztr(i)+x(ii+3)*wzt_ver
	t_ztr(i)=t_ztr(i)+x(ii+4)*wzt_time
	write(12,*)x_ztr(i),y_ztr(i),z_ztr(i),t_ztr(i)
	av1=av1+sqrt(x_ztr(i)*x_ztr(i)+y_ztr(i)*y_ztr(i))
	av2=av2+abs(z_ztr(i))
	av3=av3+abs(t_ztr(i))
end do
av1=av1/nztr
av2=av2/nztr
av3=av3/nztr
write(*,*)' source corrections:'
write(*,*)' av1=',av1,' av2=',av2,' av3=',av3
close(12)


STOP
end

