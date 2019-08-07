character*8 ar,md,line
character*1 it,itold,ppss,atest,rm,gr,ps

real xray(1000),yray(1000),zray(1000)
real resid(100),ztest(4),yprof(10,50)

integer iotper(100),iotr4(4)

real  xtop(:,:), ztop(:,:),vtop(:,:)
integer popor(:,:),obr(:,:),kod_otr(:,:,:),notr(:),ntop(:)
real ornt(20),ylevel(:)


real zper(500),votr(4),zotr(4),vaver(8)

real vertet(4,3),vtoptet(4)
real d00(1250),r00(1250)
real aa8(4,4),q8(4),a8(4)
real aa(4,4),q(4),a(4),curr(3)
integer iuz(8),nurr(8),muzel(1200)
real matr(:),matruzel(1200)
integer npar(50),nprof(10),ist1(100)


common/velmod/nnps,ypros,xpros,zpros,cfps
common/inipar/nzone,nsurf,smth,xlim1,xlim2,ylim1,ylim2,zlim1,zlim2
common/model/ar
common/itstep/itstep

allocatable xtop,ztop,vtop,matr
allocatable popor,obr,kod_otr,notr,ntop,ylevel

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

k_reject=1
rz=6371.


open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)ips		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(ps,'(i1)')ips
write(gr,'(i1)')igr

write(*,*)' execution of matr'
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
close(1)

if(ips.eq.2) att_aver=val21; add_att=val22

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
read(1,*)(ornt(i),i=1,nornt)
close(1)

orient=ornt(igr)
write(*,*)' orient=',orient
sinbase=sin(orient*per)
cosbase=cos(orient*per)

!******************************************************************
if(key_ft1_xy2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
    do i=1,10000
	    read(1,'(a8)',end=553)line
	    if(line.eq.'AREA_CEN') goto 554
    end do
    553 continue
    write(*,*)' cannot find AREA CENTER in MAJOR_PARAM.DAT!!!'
    pause
    554 read(1,*)fi0,tet0
    write(*,*)fi0,tet0
    close(1)
else
    fi0=0; tet0=0
end if

!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=533)line
	if(line.eq.'GRID_PAR') goto 534
end do
533 continue
write(*,*)' cannot find GRID_PARAMETERS in MAJOR_PARAM.DAT!!!'
pause
534 continue
read(1,*)xlim1,xlim2,dxpl
read(1,*)ylim1,ylim2,dypl
read(1,*)zlim1,zlim2,dzpl
close(1)
!******************************************************************

nsurf=0

call read_vref(ar,md)
call read_3D_mod_v(ar,md,iter)


!**************************************************************
!*    Read initial parameters of grid

nmax=0
notmax=0


open(1,file='../../../DATA/'//ar//'/'//md//'/data/lev_att'//gr//ps//'.dat')
i=0
722 i=i+1
	read(1,*,end=721)n,y
	goto 722
721 nypl=i-1
!write(*,*)' nypl=',nypl
nlev=nypl
close(1)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
do n=1,nypl
	read(1,*)nt
	if(nt.gt.nmax) nmax=nt
	do i=1,nt
		read(1,*)x,z,l
	end do
end do
close(1)

open(4,file='../../../TMP_files/tmp/otr_att'//gr//ps//'.dat')
do nur=1,nypl
	read(4,*) no
	if(no.gt.notmax) notmax=no
	read(4,*) ((ll, i=1,2), j=1,no)
end do
close(4)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr_att'//gr//ps//'.dat')
read(1,*) nobr
close(1)

nymax=nlev

ntmax=nobr

write(*,*)' nymax=',nymax,' ntmax=',ntmax
write(*,*)' nmax=',nmax,' otr nmax=',notmax


allocate(ntop(nymax),ylevel(nymax))
allocate(xtop(nmax,nymax),ztop(nmax,nymax))
allocate(vtop(nmax,nymax),popor(nmax,nymax))
allocate(kod_otr(2,notmax,nymax),notr(nymax))
allocate(obr(ntmax,2),matr(ntmax))

open(1,file='../../../DATA/'//ar//'/'//md//'/data/lev_att'//gr//ps//'.dat')
do i=1,nlev
	read(1,*)ntop(i),ylevel(i)
end do
close(1)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
do n=1,nlev
    read(1,*)nt
    ntop(n)=nt
    do i=1,ntop(n)
        read(1,*)xtop(i,n),ztop(i,n),popor(i,n)
        xx=xtop(i,n)
        yy=ylevel(n)
        zz=ztop(i,n)

        vtop(i,n)=velocity(xx,yy,zz,ips)
    end do
end do
close(1)

open(4,file='../../../TMP_files/tmp/otr_att'//gr//ps//'.dat')
do nur=1,nlev
	read(4,*) notr(nur)
	read(4,*) ((kod_otr(i,j,nur), i=1,2), j=1,notr(nur))
end do
close(4)


open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr_att'//gr//ps//'.dat')
read(1,*) nobr
read(1,*)((obr(i,j),i=1,nobr), j=1,2)
close(1)
write(*,*)' number of velocity parameters=',nobr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


nray=0
nonzer=0

!write(*,*)' nlev(1)=',nlev(1)


nr=0

open(1,file='../../../DATA/'//ar//'/'//md//'/data/atten'//ps//'.dat')
open(2,file='../../../DATA/'//ar//'/'//md//'/data/ray_paths_att'//ps//'.dat',form='binary')
open(11,file='../../../TMP_files/tmp/matr_att'//gr//ps//'.dat')

728 continue
    read(1,*,end=729)xzt,yzt,zzt,xst,yst,zst,att0,ttt
    nray=nray+1

    read(2)npr
    if(npr.eq.0)goto 728
    do i=1,npr
        read(2)xray(i),yray(i),zray(i)
    end do

    att = att0 - (add_att + att_aver*ttt)

    !write(*,*)att0,ttt,(add_att - att_aver*ttt),att

    !if(abs(att).gt.att_range) cycle

    nr=nr+1
    !if(nr.lt.71) cycle
    matr=0.
    do ipt=2,npr
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
        vmid=velocity(xmid,ymid,zmid,ips)

        xxx=xmid*cosbase+ymid*sinbase
        yyy=-xmid*sinbase+ymid*cosbase
        zzz=zmid

        do n=1,nlev-1
	        yl1=ylevel(n)
	        yl2=ylevel(n+1)
	        if ((yyy-yl1)*(yyy-yl2).le.0.) exit
        end do
        nur=n
        nur1=nur
        nur2=nur+1


        224			continue
        if((xxx-xlim1)*(xxx-xlim2).ge.0.) cycle
        if((yyy-ylim1)*(yyy-ylim2).ge.0.) cycle
        if((zzz-zlim1)*(zzz-zlim2).ge.0.) cycle


        iper=0
        do iotr=1,notr(nur1)
            np1=kod_otr(1,iotr,nur1)
            np2=kod_otr(2,iotr,nur1)
            xp1=xtop(np1,nur1)
            xp2=xtop(np2,nur1)
            if(abs(xp1-xp2).lt.0.00001) cycle
            !write(*,*)' xp1=',xp1,' xp2=',xp2,' xxx=',xxx
            if((xp1-xxx)*(xp2-xxx).gt.0.) cycle
            !write(*,*)' xp1=',xp1,' xp2=',xp2
            zp1=ztop(np1,nur1)
            zp2=ztop(np2,nur1)
            !write(*,*)' zp1=',zp1,' zp2=',zp2
            iper=iper+1
            zper(iper)=zp1+((zp2-zp1)/(xp2-xp1))*(xxx-xp1)
            iotper(iper)=iotr
            !write(*,*)iper,' zp=',zper(iper),' iotr=',iotr
            !write(*,*)
        end do

        do ip=2,iper
	        z1=zper(ip-1)
	        z2=zper(ip)
	        !write(*,*)ip,' z1=',z1,' z2=',z2,' zzz=',zzz
	        if((z1-zzz)*(z2-zzz).le.0.) goto 710
        end do
        cycle


        iper=0
        do iotr=1,notr(nur1)
	        np1=kod_otr(1,iotr,nur1)
	        np2=kod_otr(2,iotr,nur1)
	        xp1=xtop(np1,nur1)
	        xp2=xtop(np2,nur1)
	        !write(*,*)' xp1=',xp1,' xp2=',xp2
	        if(abs(xp1-xp2).lt.0.00001) cycle
	        if((xp1-xxx)*(xp2-xxx).gt.0.) cycle
	        zp1=ztop(np1,nur1)
	        zp2=ztop(np2,nur1)
	        iper=iper+1
	        zper(iper)=zp1+((zp2-zp1)/(xp2-xp1))*(xxx-xp1)
	        iotper(iper)=iotr
	        !write(*,*)' zp=',zper(iper),' iotr=',iotr
        end do


        do ip=1,iper
	        write(*,*)' zper=',zper(ip)
        end do

        write(*,*)' xxx=',xxx,' yyy=',yyy,' zzz=',zzz
        write(*,*)' zlim1=',zlim1,' zlim2=',zlim2
        write(*,*)' out of the study volume'
        pause


        710			continue
        iotr4(1)=iotper(ip-1)
        iotr4(2)=iotper(ip)
        !write(*,*)' iotr1=',iotr4(1),' iotr2=',iotr4(2)
        !write(*,*)' iper1=',iper


        iper=0
        do iotr=1,notr(nur2)
	        np1=kod_otr(1,iotr,nur2)
	        np2=kod_otr(2,iotr,nur2)
	        xp1=xtop(np1,nur2)
	        xp2=xtop(np2,nur2)
	        if(abs(xp1-xp2).lt.0.00001) cycle
	        if((xp1-xxx)*(xp2-xxx).gt.0.) cycle
	        zp1=ztop(np1,nur2)
	        zp2=ztop(np2,nur2)
	        !write(*,*)iotr,' x1=',xp1,xp2,' z1=',zp1,zp2
	        iper=iper+1
	        zper(iper)=zp1+((zp2-zp1)/(xp2-xp1))*(xxx-xp1)
	        iotper(iper)=iotr
        end do


        do ip=2,iper
	        z1=zper(ip-1)
	        z2=zper(ip)
	        if((z1-zzz)*(z2-zzz).le.0.) goto 711
        end do
        cycle

        write(*,*)' iper=',iper
        do ip=1,iper
	        write(*,*)' zper=',zper(ip)
        end do


        write(*,*)' xxx=',xxx,' yyy=',yyy,' zzz=',zzz
        write(*,*)' out of the study volume'
        pause
        711			continue
        iotr4(3)=iotper(ip-1)
        iotr4(4)=iotper(ip)
        !write(*,*)' iotr3=',iotr4(3),' iotr4=',iotr4(4)

        iinf=0
        do iotr=1,4
            nnn=nur1
            if(iotr.gt.2) nnn=nur2
            do iside=1,2
                iuzz=kod_otr(iside,iotr4(iotr),nnn)
                imatr=popor(iuzz,nnn)
                if(imatr.eq.0) cycle
                if(iinf.ne.0) then
                    do in=1,iinf
	                    if(iuz(in).eq.iuzz.and.nurr(in).eq.nnn) goto 712
                    end do
                end if
                iinf=iinf+1
                !write(*,*)iinf,iuzz,nnn,imatr
                iuz(iinf)=iuzz
                nurr(iinf)=nnn
                !write(*,*)' iuz=',iuz(1),' nurr=',nurr(1)
                vaver(iinf)=att_aver
                !write(*,*)' iuz=',iuz(1),' nurr=',nurr(1)
            712					continue
            end do
        end do

        do iotr=1,4
	        nnn=nur1
	        if(iotr.gt.2) nnn=nur2
	        ii1=kod_otr(1,iotr4(iotr),nnn)
	        x1=xtop(ii1,nnn)
	        z1=ztop(ii1,nnn)
	        v1=vtop(ii1,nnn)
	        ii2=kod_otr(2,iotr4(iotr),nnn)
	        x2=xtop(ii2,nnn)
	        z2=ztop(ii2,nnn)
	        v2=vtop(ii2,nnn)
	        !write(*,*)' v1=',v1,' v2=',v2
	        votr(iotr)=v1+((v2-v1)/(x2-x1))*(xxx-x1)
	        zotr(iotr)=z1+((z2-z1)/(x2-x1))*(xxx-x1)
        end do
        z11=zotr(1)
        z21=zotr(2)
        z12=zotr(3)
        z22=zotr(4)
        v11=votr(1)
        v21=votr(2)
        v12=votr(3)
        v22=votr(4)
        !write(*,*)' v11=',v11,' v21=',v21
        !write(*,*)' v12=',v12,' v22=',v22
        vv1=v11+((v21-v11)/(z21-z11))*(zzz-z11)
        vv2=v12+((v22-v12)/(z22-z12))*(zzz-z12)
        !write(*,*)' vv1=',vv1,' vv2=',vv2
        y1=ylevel(nur1)
        y2=ylevel(nur2)
        v000=att_aver

!			write(*,*)' v000=',v000,' iinf=',iinf

        do itop=1,iinf
            !write(*,*)' iuz=',iuz(itop),' nurr=',nurr(itop)
            !pause
            iuz0=iuz(itop)
            nur0=nurr(itop)
            imatr=popor(iuz0,nur0)
            if(imatr.eq.0) pause

            do iotr=1,4
                nnn=nur1
                if(iotr.gt.2) nnn=nur2
                !write(*,*)' iotr4(iotr)=',iotr4(iotr)
                !pause
                ii1=kod_otr(1,iotr4(iotr),nnn)
                !write(*,*)' ii1=',ii1
                x1=xtop(ii1,nnn)
                z1=ztop(ii1,nnn)
                v1=att_aver
                if(ii1.eq.iuz0.and.nnn.eq.nur0) then
	                !write(*,*)' in node:',v1,v1*1.02
	                v1=v1*1.02
                end if
                ii2=kod_otr(2,iotr4(iotr),nnn)
                x2=xtop(ii2,nnn)
                z2=ztop(ii2,nnn)
                v2=att_aver
                if(ii2.eq.iuz0.and.nnn.eq.nur0) then
	                !write(*,*)' in node:',v2,v2*1.02
	                v2=v2*1.02
                end if
                !write(*,*)' v1=',v1,' v2=',v2
                votr(iotr)=v1+((v2-v1)/(x2-x1))*(xxx-x1)
                zotr(iotr)=z1+((z2-z1)/(x2-x1))*(xxx-x1)
            end do
            !pause
            z11=zotr(1)
            z21=zotr(2)
            z12=zotr(3)
            z22=zotr(4)
            v11=votr(1)
            v21=votr(2)
            v12=votr(3)
            v22=votr(4)
            !write(*,*)' v11=',v11,' v21=',v21
            !write(*,*)' v12=',v12,' v22=',v22
            vv1=v11+((v21-v11)/(z21-z11))*(zzz-z11)
            vv2=v12+((v22-v12)/(z22-z12))*(zzz-z12)
            !write(*,*)' vv1=',vv1,' vv2=',vv2
            y1=ylevel(nur1)
            y2=ylevel(nur2)
            vfound=vv1+((vv2-vv1)/(y2-y1))*(yyy-y1)

            deltav=vfound-v000
            deltat=deltav*ds !/vmid
            !write(*,*)' deltat=',deltat
            if(vaver(itop).lt.0.00001) then
                write(*,*)' vaver=',vaver(itop)
                pause
            end if
            dmatr=deltat/(0.02*att_aver)
            matr(imatr)=matr(imatr)+dmatr
!write(*,*)' imatr=',imatr,' matr(imatr)=',matr(imatr)
            end do
        end do
        nuz=0
        do i=1,nobr
            if (abs(matr(i)).lt.1.e-10) cycle
            nuz=nuz+1
            muzel(nuz)=i
            matruzel(nuz)=matr(i)
        !write(*,*)' muzel=',muzel(nuz),' m=',matruzel(nuz)
        end do
        !pause

        if(nuz.ne.0) then
            nonzer=nonzer+nuz
            write(11,*) nuz,att
            write(11,*) (muzel(i),matruzel(i),i=1,nuz)
            if(mod(nray,50).eq.0)write(*,*)nray,nuz,att
       end if

    goto 728
    729	close(1)
close(2)

close(11)

open(11,file='../../../DATA/'//ar//'/'//md//'/data/numbers'//gr//ps//'.dat')
write(11,*)nray,nobr,nonzer
close(11)

write(*,*)' nray=',nray,' nonzer=',nonzer

stop
end
