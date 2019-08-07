character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*3 sec
character*1 it,ppss,gr,ps
real plll(600),cfff(600),zsurf(40),sumzone(40)
real sumver(600),zgv(600),verlin(600)
integer ngrd(600),kodgv(600),nbl_y(600)
integer popor(:,:),obr(:,:)

allocatable xgrd(:,:),zgrd(:,:),xnew(:),znew(:)
allocatable ival(:,:),ivzone(:),ivnew(:),popor,obr

allocatable plotray(:,:,:),plot2(:,:),sumtot(:),coeff(:),npar(:)
allocatable sumnetx(:)

one=1.d0
pi=asin(one)*2.d0
per=pi/180.d0

ndiv=2

nsurf=0

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)ips		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(ps,'(i1)')ips
write(gr,'(i1)')igr

write(*,*)' execution grid'
write(*,*)' ar=',ar,' md=',md,' ps=',ps,' gr=',gr


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
read(1,*)plotmin,plotmax	! maximal ray density, relative to average
close(1)

distmin=dzpl
zupper=zlim1
dx=dzpl/20; dz=dzpl/20


nmax_p=100000
nmax_s=100000
nvermin=5
coefmax=20


nxpl=int_best((xlim2-xlim1)/dxpl)
nypl=int_best((ylim2-ylim1)/dypl)
nzpl=int_best((zlim2-zlim1)/dzpl)


write(*,*)' nx=',nxpl,' ny=',nypl,' nz=',nzpl
write(*,*)' nmax_p=',nmax_p

allocate(plotray(nxpl,nypl,nzpl),plot2(nxpl,nzpl))
allocate(sumtot(nypl+1),coeff(nypl+1),npar(nypl+1))
allocate(sumnetx(nxpl))

allocate(popor(nmax_p,nypl+1),ival(nmax_p,nypl+1),obr(nmax_p,2))
allocate(xgrd(nmax_p,nypl+1),zgrd(nmax_p,nypl+1))
allocate(xnew(nmax_p),znew(nmax_p),ivzone(nmax_p),ivnew(nmax_p))


open(1,file='../../../COMMON/various/aver_koef.dat')
i=0
838	i=i+1
	read(1,*,end=837)plll(i),cfff(i)
!	write(*,*)plll(i),cfff(i)
	goto 838
837	close(1)
mmm=i-1
!write(*,*)' mmm=',mmm

open(1,file='../../../DATA/'//ar//'/'//md//'/data/num_att'//ps//'.dat')
read(1,*) nray
close(1)


ngrid=nmax_p
if(iiips.eq.2) ngrid=nmax_s

sumtotal=0.
nonzer=0
open(1,file='../../../TMP_files/tmp/plot_att'//gr//ps//'.dat')
do iy=1,nypl
    read(1,*)((plotray(ix,iy,iz),ix=1,nxpl),iz=1,nzpl)
    do iz=1,nzpl
        do ix=1,nxpl
            if(plotray(ix,iy,iz).lt.1.)cycle
            nonzer=nonzer+1
            sumtotal=sumtotal+plotray(ix,iy,iz)
        end do
    end do
end do
close(1)
averplot=sumtotal/nonzer

write(*,*)' aver ray lenght in one block=',averplot

do iy=1,nypl
    do iz=1,nzpl
        do ix=1,nxpl
            if(plotray(ix,iy,iz).gt.averplot*plotmax)then
                plotray(ix,iy,iz)=averplot*plotmax
            else if(plotray(ix,iy,iz).lt.averplot*plotmin)then
                plotray(ix,iy,iz)=0.
            end if
        end do
    end do
end do

summm=0
do iy=1,nypl
    sum=0
    plmax=0
    nonzer=0

    do ix=1,nxpl
        do iz=1,nzpl
            if(plotray(ix,iy,iz).lt.1.)cycle
            nonzer=nonzer+1
            plot=plotray(ix,iy,iz)
            sum=sum+plot
            if(plot.gt.plmax)plmax=plot
        end do
    end do
    nbl_y(iy)=nonzer

    sumtot(iy)=sum
    summm=summm+sum
    averplot=sumtotal/nonzer
    !write(*,*)iy,' sumtot=',sumtot(iy),' nbl_y=',nbl_y(iy)

end do
!pause

sumtot=100 * sumtot/summm
!write(*,*)' summ=',summm

do iy=1,nypl
    coeff(iy)=0
    if(sumtot(iy).lt.plll(1))cycle
    do i=2,mmm
        pl1=plll(i-1)
        pl2=plll(i)
        if((pl1-sumtot(iy))*(pl2-sumtot(iy)).gt.0.) cycle
        cf1=cfff(i-1)
        cf2=cfff(i)
        coeff(iy)=cf1+((cf2-cf1)/(pl2-pl1))*(sumtot(iy)-pl1)
        exit
    end do
end do

sum=0
do i=1,nypl
    if(sumtot(i).lt.plll(1))cycle
    sum=sum+coeff(i)
    !write(*,*)' sum=',sum
end do

isum=0
do i=1,nypl
    if(sumtot(i).eq.0)cycle
    npar(i)=coeff(i)*ngrid/sum
    isum=isum+npar(i)
    !write(*,*)' level:',i,' sum=',sumtot(i),' npar=',npar(i)
end do
!write(*,*)' Total parameters:=',isum

open(12,file='../../../TMP_files/tmp/nodes'//ps//'.dat')


do iy=1,nypl


    npar(iy)=0

    yy=ylim1+(iy-0.5)*dypl
    !call decsf(0.,yy,0.,fi0,tet0,FI,TET,h)
    write(sec,'(i3)')iy

    do ix=1,nxpl
        do iz=1,nzpl
            plot2(ix,iz)=plotray(ix,iy,iz)
        end do
    end do
    nparam=nbl_y(iy)/ndiv
    if(nparam.eq.0) cycle

    aver=sqrt(nparam*( (xlim2-xlim1) / (zlim2-zupper) ))
    nver=int(aver)
    !write(*,*)' xlim1=',xlim1,' xlim2=',xlim2
    !write(*,*)' zupper=',zupper,' zlim2=',zlim2
    !write(*,*)' nver=',nver,' nparam=',nparam
		
    kase=1
    if(abs((nver+1.)-aver).lt.abs(nver-aver))nver=nver+1
    if(nver.lt.nvermin) nver=nvermin
    sumtotal=0.
    nonzer=0

    xpl1=xlim1+0.5*dxpl
    xpl2=xlim1+(nxpl-0.5)*dxpl
    zpl1=zlim1+0.5*dzpl
    zpl2=zlim1+(nzpl-0.5)*dzpl

    averplot=sumtotal/nonzer

    do ix=1,nxpl
        sumnetx(ix)=0
        do iz=1,nzpl
            sumnetx(ix) = sumnetx(ix) + plotray(ix,iy,iz)
        end do
    end do

    do ix=1,nxpl
	    if (sumnetx(ix).ne.0.) exit
    end do
    xleft=xlim1+(ix-2)*dxpl

    do ix=nxpl,1,-1
	    if (sumnetx(ix).ne.0.) exit
    end do
    xrigt=xlim1+(ix+1)*dxpl


    dxband=dypl

    631	continue

    nver=(xrigt-xleft)/dxband
	
	
    sumtott=0.


    do iver=1,nver
        xver1=xleft+(iver-1)*dxband
        xver2=xleft+iver*dxband

        sum1=alin(xver1, nxpl,xlim1,dxpl,sumnetx)
        sum2=alin(xver2, nxpl,xlim1,dxpl,sumnetx)

        sumver(iver)=0.
        do ix=1,nxpl-1
            xn1=xlim1+(ix-0.5)*dxpl
            xn2=xlim1+(ix+0.5)*dxpl
            if(xn2.lt.xver1) cycle
            if(xn1.gt.xver2) exit
            x1=xver1
            s1=sum1
            x2=xver2
            s2=sum2
            if(xn1.gt.xver1)then
	            x1=xn1
	            s1=sumnetx(ix)
            end if
            if(xn2.lt.xver2)then
	            x2=xn2
	            s2=sumnetx(ix+1)
            end if
            dsum=(s1+s2)*(x2-x1)/2.
            sumver(iver)=sumver(iver)+dsum
        end do
        sumtott=sumtott+sumver(iver)
    end do

    ! compute the number of nodes at each band
    sumpnt=sumtott/nparam


    ngrdtot=0
    ermin=100000
    ermax=-100000
    imax=0
    imin=0
    do iver=1,nver
        ccc=sumver(iver)/sumpnt
        if(ccc.lt.1.)then
	        ngrd(iver)=0
	        cycle
        end if
        ngr=int(ccc)
        if(abs((ngr+1.)-ccc).lt.abs(ngr-ccc))ngr=ngr+1
        err=ccc-ngr
        if(err.lt.ermin)then
	        ermin=err
	        imin=iver
        end if
        if(err.gt.ermax)then
	        ermax=err
	        imax=iver
        end if
        ngrdtot=ngrdtot+ngr
        ngrd(iver)=ngr
        !		write(*,*)' ngr=',ngr
    end do
    if(imin.eq.0.and.imax.eq.0) then
	    npar(iy)=0
	    cycle
    end if

456 if(ngrdtot.lt.nparam) ngrd(imax)=ngrd(imax)+1
    if(ngrdtot.lt.nparam) ngrdtot=ngrdtot+1
    if(ngrdtot.gt.nparam) ngrd(imin)=ngrd(imin)-1
    if(ngrdtot.gt.nparam) ngrdtot=ngrdtot-1
    if(ngrdtot.eq.nparam) goto 457
    ermin=100000 
    ermax=-100000 

    do iver=1,nver 
        ccc=sumver(iver)/sumpnt 
        if(ccc.lt.1.)then
	        ngrd(iver)=0
	        cycle
        end if
        ngr=ngrd(iver) 
        err=ccc-ngr 
        if(err.lt.ermin.and.ngr.ne.0)then 
	        ermin=err 
	        imin=iver 
        end if 
        if(err.gt.ermax)then 
	        ermax=err 
	        imax=iver 
        end if 
    end do 
    goto 456

    457	continue

    ngrdtot=0
    node=0

    ! Left side of the study volume:

    node=node+1
    xgrd(node,iy)=xlim1
    zgrd(node,iy)=zlim1
    ival(node,iy)=0

    node=node+1
    xgrd(node,iy)=xlim1
    zgrd(node,iy)=zlim2
    ival(node,iy)=0

! Nodes on the vertical lines:

    do iver=1,nver
        if(ngrd(iver).eq.0) cycle
        igr=ngrd(iver)

        xver1=xleft+(iver-1)*dxband
        xver2=xleft+iver*dxband
        xnod=(xver2+xver1)/2.

        ! Top node:
        node=node+1
        xgrd(node,iy)=xnod
        zgrd(node,iy)=zlim1
        ival(node,iy)=0

        igv=0		! igv is a counter for the number of nodes in one vertical line
        zsrf1=zlim1
        zsrf2=zlim2

        sum=0
        zcur=zsrf1-dz
        123		zcur=zcur+dz
        if(zcur.gt.zsrf2) goto 36
            fun2=alin2d(xnod,zcur, xlim1,nxpl,dxpl, zlim1,nzpl,dzpl, plot2)
            dsum=dz*fun2
            sum=sum+dsum
            !write(*,*)' zcur=',zcur,' fun2=',fun2,' sum=',sum
            goto 123
        36		continue
        sumband=sum

	porog = sumband / (igr*2.)
!			write(*,*)' isurf=',isurf,' igr=',igr,' sumzone=',sumzone(isurf)
!			write(*,*)' porog=',porog
!			write(*,*)' surface1: ',zsrf1
	if(porog.lt.sumpnt/100) goto 358

	zcur=zsrf1
	nrect=0
24		sum=0.
	dist=0.
23	    zcur=zcur+dz
	    dist=dist+dz
	    if(zcur.gt.zsrf2-distmin*2) goto 358
	    fun2=alin2d(xnod,zcur, xlim1,nxpl,dxpl, zlim1,nzpl,dzpl, plot2)
	    dsum=dz*fun2
	    sum=sum+dsum
    !				write(*,*)' sum=',sum,' dsum=',dsum,' zcur=',zcur
	    if (sum.lt.porog.or.dist*2.lt.distmin) goto 23
				
    356		nrect=nrect+1
	    if(abs((nrect/2.)-int(nrect/2.)).gt.0.0000001) goto 24

	    igv=igv+1
	    zgv(igv)=zcur
	    kodgv(igv)=1
	    goto 24

    358		continue


	    do i=1,igv
		    node=node+1
		    xgrd(node,iy)=xnod
		    zgrd(node,iy)=zgv(i)

		    !if(zgv(i).gt.300..and.kodgv(i).eq.1) then
			    !write(*,*)xnod,zgv(i),kodgv(i)
			    !pause
		    !end if
		    ival(node,iy)=kodgv(i)
	    end do
	    node=node+1
	    xgrd(node,iy)=xnod
	    zgrd(node,iy)=zlim2
	    ival(node,iy)=0

    end do		


! Right side of the study volume:
    node=node+1
    xgrd(node,iy)=xlim2
    zgrd(node,iy)=zlim1
    ival(node,iy)=0


    node=node+1
    xgrd(node,iy)=xlim2
    zgrd(node,iy)=zlim2
    ival(node,iy)=0


    632	continue
    do i=2,node
	    x1=xgrd(i-1,iy)
	    x2=xgrd(i,iy)
	    if(abs(x1-x2).gt.0.000001)cycle
	    z1=zgrd(i-1,iy)
	    z2=zgrd(i,iy)
	    if(abs(z1-z2).lt.dz)then
		    if(i.eq.2.or.i.eq.node) goto 633
		    zz1=zgrd(i-2,iy)
		    zz2=zgrd(i+1,iy)
		    iv1=ival(i-2,iy)
		    iv2=ival(i+1,iy)
		    if(abs(zz1-z1).lt.coefmax.and.iv1.ne.0)ival(i-1,iy)=2
		    if(abs(zz2-z2).lt.coefmax.and.iv2.ne.0)ival(i,iy)=2
	    end if
    633		continue
    end do

    npar(iy)=node
    if(mod(iy,20).eq.0)write(*,*)' iy=',iy,' node=',node,' yy=',yy


end do
close(1)
close(12)

open(12,file='../../../DATA/'//ar//'/'//md//'/data/npar_att'//gr//ps//'.dat')
write(12,*)(npar(i),i=1,nypl)
close(12)

open(11,file='../../../TMP_files/tmp/ver_att'//gr//ps//'.dat')
open(12,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
open(13,file='../../../DATA/'//ar//'/'//md//'/data/pop_att'//gr//ps//'.dat')
open(14,file='../../../DATA/'//ar//'/'//md//'/data/obr_att'//gr//ps//'.dat')
open(15,file='../../../DATA/'//ar//'/'//md//'/data/lev_att'//gr//ps//'.dat')

npop=0
popor=0
nylevel=1

do iy=1,nypl
    nparam=npar(iy)
    if(nparam.eq.0) cycle

    nylevel=nylevel+1
		

    nnew=1
    xnew(nnew)=xgrd(1,iy)
    znew(nnew)=zgrd(1,iy)
    do i=2,nparam
        dzzz=abs(znew(nnew)-zgrd(i,iy))
        if(dzzz.gt.2..or.dzzz.lt.0.0001) then
	        nnew=nnew+1
	        xnew(nnew)=xgrd(i,iy)
	        znew(nnew)=zgrd(i,iy)
	        ivnew(nnew)=ival(i,iy)
        end if
    end do
    nparam=nnew

    iver=1
    verlin(1)=xgrd(1,iy)
    do i=2,nparam
        if(abs(xnew(i)-xnew(i-1)).gt.dx) then
	        iver=iver+1
	        verlin(iver)=xnew(i)
        end if
        if(ivnew(i).eq.0) cycle
        npop=npop+1
        !		write(*,*)npop,iy
        popor(i,nylevel)=npop
        obr(npop,1)=i
        obr(npop,2)=nylevel
    end do
    !write(*,*)nylevel,' nparam=',nparam,' npop=',npop
    inw=nparam
    yy=ylim1+(iy-0.5)*dypl

    do i=1,inw
        xx=xnew(i)
        zz=znew(i)
        if(i.ne.1)then
	        if(abs(zz-znew(i-1)).lt.0.5) zz=zz+0.1
        end if
        if(i.ne.inw)then
	        if(abs(zz-znew(i+1)).lt.0.5) zz=zz-0.1
        end if
        ivzone(i)=1
        !ivzone(i)=numzone(xx,yy,zz)
    end do



    write(11,*)iver
    write(11,*)(verlin(i),i=1,iver)
    if(nylevel.eq.2) then
        ylevel=ylim1
        write(11,*)iver
        write(11,*)(verlin(i),i=1,iver)
        write(12,*)inw
        do i=1,inw
            ipop=0
            write(12,*)xnew(i),znew(i),ipop,ivzone(i)
        end do
        write(15,*)inw,ylevel
        iii=0.
        write(13,*)(iii,i=1,inw)
    end if
    ylevel=ylim1+(iy-0.5)*dypl
    write(11,*)iver
    write(11,*)(verlin(i),i=1,iver)
    write(12,*)inw
    do i=1,inw
	    write(12,*)xnew(i),znew(i),popor(i,nylevel),ivzone(i)
    end do
    write(15,*)inw,ylevel
    !write(*,*)inw,ylevel
    write(13,*)(popor(i,nylevel),i=1,inw)

end do
!write(*,*)' npop=',npop
!pause

! The last y-level with zero values of grid
write(12,*)inw
do i=1,inw
	ipop=0
	write(12,*)xnew(i),znew(i),ipop,ivzone(i)
end do
write(11,*)iver
write(11,*)(verlin(i),i=1,iver)
iii=0.
write(13,*)(iii,i=1,inw)
ylevel=ylim2
write(15,*)inw,ylevel

write(14,*)npop
write(*,*)' number of valuable velocity parameters:',npop
write(14,*)((obr(i,j),i=1,npop), j=1,2)
!do i=1,npop
!	iuzz=obr(i,1)
!	nurrr=obr(i,2)
!	ipop=popor(iuzz,nurrr)
!	write(*,*)iuzz,nurrr,i,ipop
!	pause
!end do
close(11)
close(12)
close(13)
close(14)

stop
end