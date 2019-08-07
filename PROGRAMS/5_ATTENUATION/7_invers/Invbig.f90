character*2 itinv
character*8 ar,md,line
character*1 it,itold,ppss,rm,gr,ps

real amatruzel(3000)
integer muzel(3000)


allocatable x(:),x0(:),u(:),v(:),w(:),aaa(:),dt(:),xmod(:)
allocatable ncolrow(:),nz_row(:)

common/velmod/nnps,ypros,xpros,zpros,cfps
common/inipar/nzone,nsurf,smth,xlim1,xlim2,ylim1,ylim2,zlim1,zlim2
common/itstep/itstep


open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)ips		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(gr,'(i1)')igr
write(ps,'(i1)')ips

write(*,*)' execution of invers'
write(*,*)' ar=',ar,' md=',md,' ps=',ps,' gr=',gr

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
	if(line.eq.'ATTENUAT') goto 574
end do
573 continue
write(*,*)' cannot find ATTENUATION PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
574 continue
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)sm_vel_p, sm_vel_s
read(1,*)rg_vel_p, rg_vel_s
read(1,*)iter_lsqr
close(1)
sm_vel=sm_vel_p; if(ips.eq.2)sm_vel=sm_vel_s 
rg_vel=rg_vel_p; if(ips.eq.2)rg_vel=rg_vel_s 
!******************************************************************

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numbers'//gr//ps//'.dat')
read(1,*)nray,nvel_att,nonzer
write(*,*)nray,nvel_att,nonzer
close(1)

ncol = nvel_att
nrows=nray

open(1,file='../../../TMP_files/tmp/num_sos_att'//gr//ps//'.dat')
read(1,*)nsos_att
write(*,*)nsos_att
close(1)

nonzer = nonzer + nsos_att*2
nrows = nrows + nsos_att
nonzer=nonzer+nvel_att
nrows=nrows+nvel_att

write(*,*)' Preliminary values: N rows=',nrows,' N nonzer=',nonzer


allocate(xmod(ncol),x(ncol),x0(ncol),v(ncol),w(ncol))
allocate(nz_row(nrows),dt(nrows),u(nrows))
allocate(aaa(nonzer),ncolrow(nonzer))

open(3,file='../../../TMP_files/tmp/matr_att'//gr//ps//'.dat')

nrow=0
nonz=0
95  continue	
    read(3,*,end=96) nuz,resid
    read(3,*) (muzel(ii),amatruzel(ii),ii=1,nuz)
    nrow=nrow+1
    do ii=1,nuz
        if(muzel(ii).eq.0) then
            write(*,*)' nrow=',nrow
            write(*,*)nuz,resid
            write(*,*) (muzel(i),amatruzel(i),i=1,nuz)
            stop
        end if
    end do

    do ii=1,nuz
        mmm=muzel(ii)
        nonz=nonz+1
        aaa(nonz) = amatruzel(ii) 
        ncolrow(nonz)=mmm
    end do
    nz_row(nrow) = nuz 
    dt(nrow) = resid 

    goto 95
96 close(3)

!write(*,*)' Main matrix: nonz=',nonz,' ir=',ir


open(2,file='../../../TMP_files/tmp/sos_att'//gr//ps//'.dat',form='unformatted')
do isos=1,nsos_att
    read(2)ipar1,ipar2,dist
    !write(*,*)ipar1,ipar2,dist
    nrow=nrow+1
    dt(nrow)=0.
    nz_row(nrow)=2

    smth=sm_vel

    nonz=nonz+1
    aaa(nonz)=smth
    ncolrow(nonz) = ipar1

    nonz=nonz+1
    aaa(nonz)=-smth
    ncolrow(nonz) = ipar2

end do
close(2)

do i=1,nvel_att 
    nrow=nrow+1
    dt(nrow)=0.
    nz_row(nrow)=1
    nonz=nonz+1
    aaa(nonz)=rg_vel
    ncolrow(nonz)= i
end do

u=dt

write(*,*)' N rows=',nrow,' N columns=',ncol,' N nonzer=',nonz
write(*,*)

!*******************************************************
!*******************************************************
!*******************************************************
call pstomo(nrow,ncol,x,u,v,w,iter_lsqr,nonz,aaa,ncolrow,nz_row)
!*******************************************************
!*******************************************************
!*******************************************************

kount=0
avres=0.
avres0=0.
do irw=1,nray
    dt1=0
    do ii=1,nz_row(irw)
        kount=kount+1
        nn=ncolrow(kount)
        dt1=dt1+aaa(kount)*x(nn)
        !write(*,*)nn,kount,aaa(kount),x(nn)
    end do
    !write(*,*)irw,dt(irw),dt1
    !pause
    avres=avres+abs(dt(irw)-dt1)
    avres0=avres0+abs(dt(irw))
end do
avres=avres/nray
avres0=avres0/nray
reduct=((avres0-avres)/avres0)*100.

write(*,*)'_____________________________________________________________'
write(*,*)' avres0=',avres0,' avres=',avres,' red=',reduct
write(*,*)'_____________________________________________________________'

open(11,file='../../../DATA/'//ar//'/'//md//'/data/att_result_'//gr//ps//'.dat')
do i=1,nvel_att
    write(11,*)x(i)
end do
close(11)


STOP
end

