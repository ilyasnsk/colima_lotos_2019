character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*3 sec
character*1 it,ppss,gr
real plll(600),cfff(600),zsurf(40),sumzone(40)
real sumver(600),zgv(600),verlin(600)
integer ngrd(600),kodgv(600),nbl_y(600)
integer popor(:,:),obr(:,:),nrps(2)

allocatable xgrd(:,:),zgrd(:,:),xnew(:),znew(:)
allocatable ival(:,:),ivzone(:),ivnew(:),obr

allocatable plotray(:,:,:),plot2(:,:),sumtot(:),coeff(:),npar(:)
allocatable sumnetx(:),popor,kod_xy(:,:,:),kall_xy(:,:,:)

one=1.d0
pi=asin(one)*2.d0
per=pi/180.d0

ndiv=2

nsurf=0

open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)iter		! code of the grid
read(1,*)igr		! code of the grid
close(1)
write(it,'(i1)')iter
write(gr,'(i1)')igr

write(*,*)
write(*,*)' *************************************************'
write(*,*)' Construction of grid for parameterization'
write(*,*)' ar=',ar,' md=',md,' it=',it,' gr=',gr

call read_topo(ar)

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
dzmini=dzpl/20

nxpl=int_best((xlim2-xlim1)/dxpl)
nypl=int_best((ylim2-ylim1)/dypl)
nzpl=int_best((zlim2-zlim1)/dzpl)


write(*,*)' nx=',nxpl,' ny=',nypl,' nz=',nzpl

allocate(plotray(nxpl,nypl,nzpl),plot2(nxpl,nzpl))
allocate(sumtot(nypl+1),coeff(nypl+1),npar(nypl+1))
allocate(sumnetx(nxpl),kod_xy(nxpl,nypl,2),kall_xy(nxpl,nypl,2))

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray'//it//'.dat')
read(1,*) nrps(1),nrps(2)
write(*,*) nrps(1),nrps(2)
close(1)

do iiips=1,2
    write(*,*)' iiips=',iiips
    if(nrps(iiips).eq.0) cycle
    kod_xy=0
    kall_xy=0
    npar=0

    write(ppss,'(i1)')iiips

    ngrid=nmax_p
    if(iiips.eq.2) ngrid=nmax_s

    sumtotal=0.
    nonzer=0
    open(1,file='../../../TMP_files/tmp/plotray'//ppss//gr//'.dat',form='binary')
    do iy=1,nypl
        read(1)((plotray(ix,iy,iz),ix=1,nxpl),iz=1,nzpl)
        !write(*,'(10f7.2)')((plotray(ix,iy,iz),ix=1,nxpl),iz=1,nzpl)
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

                if(plotray(ix,iy,iz).lt.1.)cycle
                    !write(*,*)plotray(ix,iy,iz),averplot*plotmin
                if(plotray(ix,iy,iz).gt.averplot*plotmax)then
                    plotray(ix,iy,iz)=averplot*plotmax
                else if(plotray(ix,iy,iz).lt.averplot*plotmin)then
                    plotray(ix,iy,iz)=0.
                end if
                !write(*,*)ix,iy,iz,plotray(ix,iy,iz)
           end do
        end do
    end do
    !!! Rays per node:
    rays_per_node = averplot
    !write(*,*)' averplot=',averplot,' nypl=',nypl
    !write(*,*)' plotray(50,62,15)=',plotray(50,62,15)
    summm=0
    nonz_tot=0
    nonz_max=0
    do iy=1,nypl
        sum=0
        plmax=0
        nonzer=0
        !write(*,*)averplot/1000.
        !write(*,'(10f7.2)')((plotray(ix,iy,iz),ix=1,nxpl),iz=1,nzpl)
       ! pause

        do ix=1,nxpl
            do iz=1,nzpl
                !if(plotray(ix,iy,iz).lt.0.000001)cycle
                !write(*,*)ix,iy,iz,plotray(ix,iy,iz)
                if(plotray(ix,iy,iz) .lt. averplot/1000. )cycle
                nonzer=nonzer+1
                nonz_tot=nonz_tot+1
                plot=plotray(ix,iy,iz)
                sum=sum+plot
                if(plot.gt.plmax)plmax=plot
            end do
        end do
        nbl_y(iy)=nonzer
        if (nonzer.gt.nonz_max) nonz_max=nonzer

        sumtot(iy)=sum
        summm=summm+sum
        !write(*,*)iy,' sumtot=',sumtot(iy),' nbl_y=',nbl_y(iy)

    end do
    !pause
    nonz_max=nonz_max*2
    write(*,*)' nonz_tot=',nonz_tot,' nonz_max=',nonz_max

    allocate(popor(nonz_max,nypl+1),ival(nonz_max,nypl+1),obr(nonz_tot,2))
    allocate(xgrd(nonz_max,nypl+1),zgrd(nonz_max,nypl+1)) 
    allocate(xnew(nonz_tot),znew(nonz_tot))


    npar=0

    xgrd(1,1)=xlim1+0.5*dxpl; zgrd(1,1)=zlim1; ival(1,1)=0
    xgrd(2,1)=xlim1+0.5*dxpl; zgrd(2,1)=zlim2; ival(2,1)=0

    xgrd(3,1)=xlim2-0.5*dxpl; zgrd(3,1)=zlim1; ival(3,1)=0
    xgrd(4,1)=xlim2-0.5*dxpl; zgrd(4,1)=zlim2; ival(4,1)=0

    npar(1)=4

    do iy=2,nypl-1

        nod_plane=int_best( sumtot(iy) / rays_per_node)
        porog = rays_per_node/2
        yy=ylim1+(iy-0.5)*dypl
        !write(*,*)' iy=',iy,' yy=',yy,' nod_plane=',nod_plane
        node = 0

        if(nod_plane.lt.1) cycle


        node=node+1
        xgrd(node,iy)=xlim1+0.5*dxpl
        zgrd(node,iy)=zlim1
        ival(node,iy)=0

        node=node+1
        xgrd(node,iy)=xlim1+0.5*dxpl
        zgrd(node,iy)=zlim2
        ival(node,iy)=0
    
        n_plane=0

        do ix=1,nxpl
            xx=xlim1+(ix-0.5)*dxpl

            !xx=0; yy=0
            call decsf(xx,yy,0.,fi0,tet0,ff,tt,h)
            
            ztopo=h_lim(ff,tt)
            !write(*,*)' ff=',ff,' tt=',tt,' ztopo=',ztopo
            !stop

            plot_on_ver=0
            do iz=1,nzpl
                plot_on_ver = plot_on_ver + plotray(ix,iy,iz)
            end do
            node_on_ver = int_best (plot_on_ver / rays_per_node)
            n_ver=0
           !write(*,*)' plot_on_ver=',plot_on_ver,' node_on_ver=',node_on_ver

            if (node_on_ver.lt.1) cycle

            !write(*,*)' xx=',xx,' yy=',yy,' ztopo=',ztopo

            node=node+1
            xgrd(node,iy)=xx
            zgrd(node,iy)=zlim1
            ival(node,iy)=0

            zcur=zlim1
            dcur=0
            ncur=0
            dzcur=0
            ngood=0

    91      zcur=zcur+dzmini
            dzcur=dzcur+dzmini
            if(zcur.ge.zlim2) goto 94

            if(zcur.lt.ztopo) goto 91

            do iz=1,nzpl
                z1=zlim1+(iz-1)*dzpl
                z2=zlim1+iz*dzpl
                if((zcur-z1)*(zcur-z2).le.0.) goto 93
            end do
            goto 94

    93      continue
            dcur = dcur + plotray(ix,iy,iz) * (dzmini/dzpl)

            !write(*,*)' zcur=',zcur,' dzcur=',dzcur,' dcur=',dcur,' porog=',porog
            if(dcur.lt.porog) goto 91
            if(dzcur.lt.distmin) goto 91
            ncur=ncur+1
            if(mod(ncur,2).ne.0) then
                node=node+1
                if(ngood.eq.0) then
                    kod_xy(ix,iy,1)=node
                    ngood=1
                end if
                n_ver=n_ver+1
                xgrd(node,iy)=xx
                zgrd(node,iy)=zcur
                ival(node,iy)=1
                !write(*,*)node,xgrd(node,iy),zgrd(node,iy),dzcur
                dzcur=0
            end if
            dcur=0
            goto 91

    94      continue
            if(ngood.eq.1) then
                kod_xy(ix,iy,2)=node
                !write(*,*)' ix=',ix,' iy=',iy,' kod1=',kod_xy(ix,iy,1),' kod2=',kod_xy(ix,iy,2)
                !pause
            end if


            if(n_ver.eq.0) then
                node=node-1
            else
                node=node+1
                n_plane=n_plane+1
                xgrd(node,iy)=xx
                zgrd(node,iy)=zlim2
                ival(node,iy)=0
            end if

        end do

        if(n_plane.eq.0) then 
            node = node - 2
        else
            node=node+1
            xgrd(node,iy)=xlim2-0.5*dxpl
            zgrd(node,iy)=zlim1
            ival(node,iy)=0


            node=node+1
            xgrd(node,iy)=xlim2-0.5*dxpl
            zgrd(node,iy)=zlim2
            ival(node,iy)=0
        end if

        npar(iy)=node
        if(mod(iy,5).eq.0)write(*,*)' iy=',iy,' node=',node,' yy=',yy
    end do

    node=0
    node=node+1
    xgrd(node,nypl)=xlim1+0.5*dxpl; zgrd(node,nypl)=zlim1; ival(node,nypl)=0
    node=node+1
    xgrd(node,nypl)=xlim1+0.5*dxpl; zgrd(node,nypl)=zlim2; ival(node,nypl)=0

    node=node+1
    xgrd(node,nypl)=xlim2-0.5*dxpl; zgrd(node,nypl)=zlim1; ival(node,nypl)=0
    node=node+1
    xgrd(node,nypl)=xlim2-0.5*dxpl; zgrd(node,nypl)=zlim2; ival(node,nypl)=0

    npar(nypl)=4


    open(12,file='../../../DATA/'//ar//'/'//md//'/data/npar'//ppss//gr//'.dat')
    write(12,*)(npar(i),i=1,nypl)
    close(12)

    open(12,file='../../../DATA/'//ar//'/'//md//'/data/kod_xy'//ppss//gr//'.dat')
    write(12,*)nxpl,nypl
    do iy=1,nypl
        write(12,*)((kod_xy(ix,iy,i2),i2=1,2),ix=1,nxpl)
    end do
    close(12)


    open(12,file='../../../DATA/'//ar//'/'//md//'/data/gr'//ppss//gr//'.dat')
    open(13,file='../../../DATA/'//ar//'/'//md//'/data/pop'//ppss//gr//'.dat')
    open(15,file='../../../DATA/'//ar//'/'//md//'/data/levinfo'//ppss//gr//'.dat')

    npop=0
    popor=0

    nylevel=0

    do iy=1,nypl

        yy=ylim1+(iy-0.5)*dypl
        write(12,*)npar(iy),yy
        write(15,*)npar(iy),yy

        if(npar(iy).eq.0) cycle

        nylevel = nylevel+1
        kall_xy(1,iy,1)=1
        xold=xlim1+0.5*dxpl
        ixold=1
        do inod=1,npar(iy)
            xxx=xgrd(inod,iy)
            !write(*,*)' xxx=',xxx
            if(abs(xxx-xold).gt.dxpl/100.) then
                do ix=1,nxpl
                    x=xlim1+(ix-0.5)*dxpl
                    if(abs(xxx-x).lt.dxpl/100.) exit
                end do
                kall_xy(ixold,iy,2)=inod-1
                kall_xy(ix,iy,1)=inod
                ixold=ix
                xold=xxx
            end if
	    if(ival(inod,iy).eq.0) cycle
            npop=npop+1
            popor(inod,iy)=npop
            obr(npop,1)=inod
            obr(npop,2)=iy
        end do
        kall_xy(nxpl,iy,2)=npar(iy)

        !write(*,*)((kall_xy(ix,iy,i2),i2=1,2),ix=1,nxpl)

        do inod=1,npar(iy)
	        write(12,*)xgrd(inod,iy),zgrd(inod,iy),popor(inod,iy)
        end do

        write(13,*)(popor(i,iy),i=1,npar(iy))
    end do
    close(15)
    close(12)
    close(13)


    open(14,file='../../../DATA/'//ar//'/'//md//'/data/obr'//ppss//gr//'.dat')
    write(14,*)npop
    write(*,*)' number of valuable velocity parameters:',npop
    write(14,*)((obr(i,j),i=1,npop), j=1,2)



    open(12,file='../../../DATA/'//ar//'/'//md//'/data/kall_xy'//ppss//gr//'.dat')
    write(12,*)nxpl,nypl
    do iy=1,nypl
        write(12,*)((kall_xy(ix,iy,i2),i2=1,2),ix=1,nxpl)
    end do
    close(12)


    deallocate(popor,ival,obr)
    deallocate(xgrd,zgrd) 
    deallocate(xnew,znew)


end do



stop
end