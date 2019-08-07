character*8 ar,md,line
character*1 it,itold,ppss,atest,rm,gr

real xray(1000),yray(1000),zray(1000)
real resid(100),ztest(4),yprof(10,50)

integer iotper(100),iotr4(4),nlev(2),nobr(2),nparam(2),nrps(2)

integer popor(:,:,:),obr(:,:,:),kod_otr(:,:,:,:),notr(:,:),ntop(:,:)
real ornt(20),ylevel(:,:)


real zper(500),votr(4),zotr(4),vaver(8)

real vertet(4,3),vtoptet(4)
real d00(1250),r00(1250)
real aa8(4,4),q8(4),a8(4)
real aa(4,4),q(4),a(4),curr(3)
integer iuz(8),nurr(8),muzel(1200)
real amatruzel(1200)
integer npar(50),nprof(10),ist1(100)
real xnod8(8),ynod8(8),znod8(8),dvnod8(8)
integer knod8(8)


common/velmod/nnps,ypros,xpros,zpros,cfps
common/inipar/nzone,nsurf,smth,xlim1,xlim2,ylim1,ylim2,zlim1,zlim2
common/model/ar
common/itstep/itstep

allocatable xtop(:,:,:), ztop(:,:,:),vtop(:,:,:),amatr(:),kall_xy(:,:,:,:)
allocatable popor,obr,kod_otr,notr,ntop,ylevel

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0

k_reject=1
rz=6371.


open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)iter		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(it,'(i1)')iter
write(gr,'(i1)')igr

write(*,*)
write(*,*)' ***********************************************'
write(*,*)' Perform MATRIX calculation: '
write(*,*)' ar=',ar,' md=',md,' it=',it,' gr=',gr

call read_vref(ar,md)
call read_3D_mod_v(ar,md,iter-1)


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
xlimm1=xlim1+0.5*dxpl
xlimm2=xlim2-0.5*dxpl
ylimm1=ylim1+0.5*dypl
ylimm2=ylim2-0.5*dypl

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray'//it//'.dat')
read(1,*) nrps(1),nrps(2)
!write(*,*) nrps(1),nrps(2)
close(1)

nmax=0
do iiips=1,2
    if(nrps(iiips).eq.0) cycle

    write(ppss,'(i1)')iiips

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/kall_xy'//ppss//gr//'.dat')
    read(1,*)nxpl,nypl
    write(*,*)nxpl,nypl
    if(iiips.eq.1) allocate(kall_xy(nxpl,nypl,2,2))
    do iy=1,nypl
        read(1,*)((kall_xy(ix,iy,i2,iiips),i2=1,2),ix=1,nxpl)
    end do
    close(1)

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr'//ppss//gr//'.dat')
    read(1,*) nparam(iiips)
    close(1)

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr'//ppss//gr//'.dat')
    do n=1,nypl
        read(1,*)nt
        !write(*,*)n,nt
        if(nt.gt.nmax) nmax=nt
        if(nt.eq.0) cycle
        do i=1,nt
	        read(1,*)x,z,l
        end do
    end do
    close(1)

end do

!write(*,'(10i5)')(kall_xy(ix,43,1,2),kall_xy(ix,43,2,2),ix=1,nxpl)

nparmax=nparam(1)
if(nparam(2).gt.nparmax) nparmax=nparam(2)
!write(*,*)' nmax=',nmax,' nparmax=',nparmax,' nparam=',nparam(1),nparam(2)

allocate(ntop(nypl,2),ylevel(nypl,2))
allocate(xtop(nmax,nypl,2),ztop(nmax,nypl,2))
allocate(vtop(nmax,nypl,2),popor(nmax,nypl,2))
allocate(obr(nparmax,2,2),amatr(nparmax))


obr=0
do iiips=1,2
    if(nrps(iiips).eq.0) cycle
    write(ppss,'(i1)')iiips

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr'//ppss//gr//'.dat')
    nparam(iiips)=0
    do iy=1,nypl
        read(1,*)nt,yy
        !write(*,*)nt,yy
        ntop(iy,iiips)=nt
        ylevel(iy,iiips)=yy
        if(nt.eq.0) cycle
        do inod=1,ntop(iy,iiips)
            read(1,*)xtop(inod,iy,iiips),ztop(inod,iy,iiips),popor(inod,iy,iiips)
            if(popor(inod,iy,iiips).ne.0) then
                nparam(iiips)=nparam(iiips)+1
                obr(nparam(iiips),1,iiips)=inod
                obr(nparam(iiips),2,iiips)=iy
            end if
            xx=xtop(inod,iy,iiips)
            yy=ylevel(iy,iiips)
            zz=ztop(inod,iy,iiips)
            vtop(inod,iy,iiips)=velocity(xx,yy,zz,iiips)
        end do
    end do
    close(1)
    write(*,*)' number of velocity parameters=',nparam(iiips)
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(11,file='../../../TMP_files/tmp/matr'//it//gr//'.dat',form='binary')

nzt=0
nray=0
nrp=0
nrs=0
nrr=0
nonzer=0
nonz_p=0
nonz_s=0

!write(*,*)' nlev(1)=',nlev(1)


nr=0
open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays'//it//'.dat')
open(2,file='../../../TMP_files/tmp/ray_paths_'//it//'.dat',form='binary')

728 continue
    if(key_true1.eq.1) read(1,end=729)xtrue,ytrue,ztrue
    read(1,*,end=729)xzt,yzt,depth,nkrat
    !write(*,*)xzt,yzt,depth,nkrat
    nzt=nzt+1
    zzt = flat_sph(key_flat1,xzt,yzt,depth)

    vzt1=velocity(xzt,yzt,zzt,1)
    vzt2=velocity(xzt,yzt,zzt,2)

    if(nkrat.eq.0) goto 728

    ik=0
    do ikr=1,nkrat
        read(1,*,end=729)ips,ist,tobs,tmod
        read(2)npr
        !write(*,*)npr
        if(npr.eq.0)cycle
        do i=1,npr
	        read(2)xray(i),yray(i),zray(i)
        end do

        nr=nr+1
!if(nr.ne.31407) cycle
!do i=1,npr
!	write(*,*)xray(i),yray(i),zray(i)
!end do
!write(*,*)nr,ist,ips,tobs,tmod



        s0=1/vzt1
        if(ips.ne.1) s0=1/vzt2
        x1=xray(1)
        y1=yray(1)
        z1=zray(1)
        x2=xray(2)
        y2=yray(2)
        z2=zray(2)

        hor=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
        tot=sqrt(hor*hor+(z2-z1)*(z2-z1))


        !write(*,*)' hor=',hor,' tot=',tot
        if(hor.gt.0.000001) then
            cosbe=(x2-x1)/hor
            sinbe=(y2-y1)/hor
            !write(*,*)' dx=',x2-x1,' dy=',y2-y1,' dz=',z2-z1
            cosal=(z2-z1)/tot
            sinal=hor/tot
            dtdx=cosbe*sinal*s0
            dtdy=sinbe*sinal*s0
            dtdz=cosal*s0
        else
            dtdx=0
            dtdy=0
            dtdz=s0
        end if

!write(*,*)' dtdx=',dtdx,' dtdy=',dtdy,' dtdz=',dtdz

        res=tobs-tmod

        amatr=0.
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

!            ztopo=topo_xy(xmid,ymid)
!            if (zmid.lt.ztopo) cycle

            vvv0=velocity(xmid,ymid,zmid,ips)

            xxx=xmid*cosbase+ymid*sinbase
            yyy=-xmid*sinbase+ymid*cosbase
            zzz=zmid

            if((xxx-xlimm1)*(xxx-xlimm2).ge.0.) cycle
            if((yyy-ylimm1)*(yyy-ylimm2).ge.0.) cycle
            if((zzz-zlim1)*(zzz-zlim2).ge.0.) cycle

    !        yyy=0.2; xxx=0.2; zzz=4
    !        vvv0=velocity(xxx,yyy,zzz)

    !        write(*,*)' xxx=',xxx,' yyy=',yyy,' zzz=',zzz,' vvv0=',vvv0

            !if(ipt.lt.41) cycle
            !write(*,'(i5,3f8.3)')ipt,xxx,yyy,zzz
            do iy=1,nypl-1
	            yl1=ylevel(iy,ips)
	            yl2=ylevel(iy+1,ips)
	            if ((yyy-yl1)*(yyy-yl2).le.0.) goto 335
            end do
            cycle
    335     nur=iy
            nur1=iy
            nur2=iy+1

    35      continue
            if(ntop(nur1,ips).eq.0) then
                nur1=nur1-1
                goto 35
            end if

    36      continue
            if(ntop(nur2,ips).eq.0) then
                nur2=nur2+1
                goto 36
            end if


            do ix=1,nxpl-1
                x1=xlim1+(ix-0.5)*dxpl
                x2=xlim1+(ix+0.5)*dxpl
                !write(*,*)ix,x1,x2,xxx
	        if ((xxx-x1)*(xxx-x2).le.0.) goto 336
            end do
            cycle
    336     continue

            ix11=ix
    37      continue
if(ix11.eq.0.or.ix11.gt.nxpl.or.nur1.eq.0.or.nur1.gt.nypl) then
    write(*,*)' PROBLEM:',ix,x1,x2,xxx
    write(*,*)' ix11=',ix11,' ix=',ix
    write(*,*)' ix11=',ix11,' nur1=',nur1,' xxx=',xxx,' yyy=',yyy
    stop
end if
            if(kall_xy(ix11,nur1,1,ips).eq.0) then
                ix11=ix11-1
                goto 37
            end if

            ix21=ix+1
    38      continue
if(ix21.eq.0.or.ix21.gt.nxpl.or.nur1.eq.0.or.nur1.gt.nypl) then
    write(*,*)' ix21=',ix21,' nur1=',nur1,' xxx=',xxx,' yyy=',yyy
    cycle
end if
            if(kall_xy(ix21,nur1,1,ips).eq.0) then
                ix21=ix21+1
                goto 38
            end if

            ix12=ix
    39      continue
if(ix12.eq.0.or.ix12.gt.nxpl.or.nur2.eq.0.or.nur2.gt.nypl) then
    write(*,*)' ix12=',ix12,' nur2=',nur2,' xxx=',xxx,' yyy=',yyy
    stop
end if
            if(kall_xy(ix12,nur2,1,ips).eq.0) then
                ix12=ix12-1
                goto 39
            end if

            ix22=ix+1
    40      continue
if(ix22.eq.0.or.ix22.gt.nxpl.or.nur2.eq.0.or.nur2.gt.nypl) then
    write(*,*)' nxpl=',nxpl,' nypl=',nypl
    write(*,*)' ix22=',ix22,' nur2=',nur2,' xxx=',xxx,' yyy=',yyy
    stop
end if
            if(kall_xy(ix22,nur2,1,ips).eq.0) then
                ix22=ix22+1
                goto 40
            end if

    !        write(*,*)' nur1=',nur1,' ix11=',ix11,' ix21=',ix21
    !        write(*,*)' nur2=',nur2,' ix12=',ix12,' ix22=',ix22

    if(ix11.eq.0.or.ix11.gt.nxpl.or.nur1.eq.0.or.nur1.gt.nypl) then
        write(*,*)' ix11=',ix11,' nur1=',nur1,' xxx=',xxx,' yyy=',yyy
        stop
    end if

            nv11_1=kall_xy(ix11,nur1,1,ips)
            nv11_2=kall_xy(ix11,nur1,2,ips)
            do iv11=nv11_1,nv11_2-1
                x111=xtop(iv11,nur1,ips)
                z111=ztop(iv11,nur1,ips)
                x112=xtop(iv11+1,nur1,ips)
                z112=ztop(iv11+1,nur1,ips)
                if((zzz-z111)*(zzz-z112).le.0.) exit
            end do
            xnod8(1)=x111; ynod8(1)=ylevel(nur1,ips); znod8(1)=z111; knod8(1)=popor(iv11,nur1,ips)
            xnod8(2)=x112; ynod8(2)=ylevel(nur1,ips); znod8(2)=z112; knod8(2)=popor(iv11+1,nur1,ips)

    if(ix21.eq.0.or.ix21.gt.nxpl.or.nur1.eq.0.or.nur1.gt.nypl) then
        write(*,*)' ix21=',ix21,' nur1=',nur1,' xxx=',xxx,' yyy=',yyy
        stop
    end if
            nv21_1=kall_xy(ix21,nur1,1,ips)
            nv21_2=kall_xy(ix21,nur1,2,ips)
            do iv21=nv21_1,nv21_2-1
                x211=xtop(iv21,nur1,ips)
                z211=ztop(iv21,nur1,ips)
                x212=xtop(iv21+1,nur1,ips)
                z212=ztop(iv21+1,nur1,ips)
                if((zzz-z211)*(zzz-z212).le.0.) exit
            end do
            xnod8(3)=x211; ynod8(3)=ylevel(nur1,ips); znod8(3)=z211; knod8(3)=popor(iv21,nur1,ips)
            xnod8(4)=x212; ynod8(4)=ylevel(nur1,ips); znod8(4)=z212; knod8(4)=popor(iv21+1,nur1,ips)

    if(ix12.eq.0.or.ix12.gt.nxpl.or.nur2.eq.0.or.nur2.gt.nypl) then
        write(*,*)' ix12=',ix12,' nur2=',nur2,' xxx=',xxx,' yyy=',yyy
        stop
    end if
            nv12_1=kall_xy(ix12,nur2,1,ips)
            nv12_2=kall_xy(ix12,nur2,2,ips)
            do iv12=nv12_1,nv12_2-1
                x121=xtop(iv12,nur2,ips)
                z121=ztop(iv12,nur2,ips)
                x122=xtop(iv12+1,nur2,ips)
                z122=ztop(iv12+1,nur2,ips)
                if((zzz-z121)*(zzz-z122).le.0.) exit
            end do
            xnod8(5)=x121; ynod8(5)=ylevel(nur2,ips); znod8(5)=z121; knod8(5)=popor(iv12,nur2,ips)
            xnod8(6)=x122; ynod8(6)=ylevel(nur2,ips); znod8(6)=z122; knod8(6)=popor(iv12+1,nur2,ips)

    if(ix22.eq.0.or.ix22.gt.nxpl.or.nur2.eq.0.or.nur2.gt.nypl) then
        write(*,*)' nxpl=',nxpl,' nypl=',nypl
        write(*,*)' ix22=',ix22,' nur2=',nur2,' xxx=',xxx,' yyy=',yyy
        stop
    end if

            nv22_1=kall_xy(ix22,nur2,1,ips)
            nv22_2=kall_xy(ix22,nur2,2,ips)
            do iv22=nv22_1,nv22_2-1
                x221=xtop(iv22,nur2,ips)
                z221=ztop(iv22,nur2,ips)
                x222=xtop(iv22+1,nur2,ips)
                z222=ztop(iv22+1,nur2,ips)
                if((zzz-z221)*(zzz-z222).le.0.) exit
            end do
            xnod8(7)=x221; ynod8(7)=ylevel(nur2,ips); znod8(7)=z221; knod8(7)=popor(iv22,nur2,ips)
            xnod8(8)=x222; ynod8(8)=ylevel(nur2,ips); znod8(8)=z222; knod8(8)=popor(iv22+1,nur2,ips)

    !        do i8=1,8
    !            write(*,*)xnod8(i8),ynod8(i8),znod8(i8),knod8(i8)
    !        end do
            do i8=1,8
                if(knod8(i8).eq.0)cycle
                dvnod8=0
                dvnod8(i8)=1
                dv11=dvnod8(1)+((dvnod8(2)-dvnod8(1))/(znod8(2)-znod8(1))) * (zzz-znod8(1))
                dv21=dvnod8(3)+((dvnod8(4)-dvnod8(3))/(znod8(4)-znod8(3))) * (zzz-znod8(3))
                dv12=dvnod8(5)+((dvnod8(6)-dvnod8(5))/(znod8(6)-znod8(5))) * (zzz-znod8(5))
                dv22=dvnod8(7)+((dvnod8(8)-dvnod8(7))/(znod8(8)-znod8(7))) * (zzz-znod8(7))
                !write(*,*)dv11,dv21,dv12,dv22

                if(abs(znod8(2)-znod8(1)).lt.0.00001)then
                    write(*,*)' znod8(1)=',znod8(1),' znod8(2)=',znod8(2)
                    stop
                end if
                if(abs(znod8(4)-znod8(3)).lt.0.00001)then
                    write(*,*)' znod8(3)=',znod8(3),' znod8(4)=',znod8(4)
                    stop
                end if
                if(abs(znod8(6)-znod8(5)).lt.0.00001)then
                    write(*,*)' znod8(5)=',znod8(5),' znod8(6)=',znod8(6)
                    stop
                end if
                if(abs(znod8(8)-znod8(7)).lt.0.00001)then
                    write(*,*)' znod8(7)=',znod8(7),' znod8(8)=',znod8(8)
                    stop
                end if

                dvv1=dv11+((dv21-dv11)/(xnod8(3)-xnod8(1))) * (xxx-xnod8(1))
                dvv2=dv12+((dv22-dv12)/(xnod8(7)-xnod8(5))) * (xxx-xnod8(5))
                !write(*,*)dvv1,dvv2

                if(abs(xnod8(3)-xnod8(1)).lt.0.00001)then
                    write(*,*)' xxx=',xxx,' yyy=',yyy,' zzz=',zzz
                    write(*,*)' ix11=',ix11,' ix21=',ix21,' nur1=',nur1,' ips=',ips
                    write(*,*)' kall_xy(ix11,nur1,1,ips)',kall_xy(ix11,nur1,1,ips),' kall_xy(ix21,nur1,1,ips)',kall_xy(ix21,nur1,1,ips)
                    write(*,*)' nv11_1=',nv11_1,' nv21_1=',nv21_1
                    write(*,*)' knod8(1)=',knod8(1),' knod8(3)=',knod8(3)
                    write(*,*)' xnod8(1)=',xnod8(1),' xnod8(3)=',xnod8(3)
                    write(*,*)' ynod8(1)=',ynod8(1),' xnod8(3)=',ynod8(3)
                    write(*,*)' znod8(1)=',znod8(1),' znod8(3)=',znod8(3)
                    stop
                end if
                if(abs(xnod8(7)-xnod8(5)).lt.0.00001)then
                    write(*,*)' xnod8(5)=',xnod8(5),' xnod8(7)=',xnod8(7)
                    stop
                end if

                dvvv=dvv1+((dvv2-dvv1)/(ynod8(5)-ynod8(1))) * (yyy-ynod8(1))

                dtdv=-dvvv*ds/(vvv0*vvv0)
                amatr(knod8(i8))=amatr(knod8(i8))+dtdv
                !write(*,*)' i8=',i8,' dvvv=',dvvv,' dtdv=',dtdv
                !write(*,*)
            end do
        end do

        nuz=0
        do i=1,nparam(ips)
            if (abs(amatr(i)).lt.1.e-10) cycle
            nuz=nuz+1
            muzel(nuz)=i
            amatruzel(nuz)=amatr(i)
            !write(*,*)' muzel=',muzel(nuz),' m=',amatruzel(nuz)
        end do

        nray=nray+1
        nonzer=nonzer+nuz
        if(ips.eq.1) then
            nrp=nrp+1
            nonz_p=nonz_p+nuz
        else
            nrs=nrs+1
            nonz_s=nonz_s+nuz
        end if

        if(mod(nray,1000).eq.0)write(*,'(4i6,f8.3)')nray,nzt,ips,nuz,res
        !write(*,*) nuz,res,ist,ips,nzt
        write(11) nuz,res,ist,ips,nzt
        write(11)dtdx,dtdy,dtdz
        if(nuz.ne.0) then
            write(11) (muzel(i),amatruzel(i),i=1,nuz)
!do i=1,nuz
!   write(*,*) muzel(i),matruzel(i)
!end do
        end if
    end do

	goto 728
729	close(1)
close(2)

close(11)

open(11,file='../../../DATA/'//ar//'/'//md//'/data/numbers'//it//gr//'.dat')
write(11,*)nray,nparam(1),nparam(2),nzt,nonzer
write(11,*)nrp,nonz_p
write(11,*)nrs,nonz_s
close(11)

write(*,*)' Number of rays=',nray,' Number of events=',nzt
write(*,*)' Numbers of P and S parameters=',nparam(1),nparam(2),' Nonzero elements in the matrix=',nonzer
write(*,*)' For P-data: nrp=',nrp,' nonz_p=',nonz_p
write(*,*)' For S-data: nrs=',nrs,' nonz_s=',nonz_s

stop
end
