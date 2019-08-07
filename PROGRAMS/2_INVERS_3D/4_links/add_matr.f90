character*8 ar,md
character*1 it,ppss,gr
integer ntop(200)
real  xtop(:,:), ztop(:,:),ylevel(200),zver1(100),zver2(100)
integer popor(:,:),nray(2)

integer links(2,1000000)

allocatable xtop,ztop,popor, kod_xy(:,:,:)


ioddeven=0

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)itstep
read(1,*)igr
close(1)

write(*,*)
write(*,*)' ******************************************************'
write(*,*)' Define LINKS between nodes'
write(*,*)' ar=',ar,' md=',md,' gr=',gr


itdvmod=itstep
!itstep=3
write(it,'(i1)')itstep
write(gr,'(i1)')igr

open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray'//it//'.dat')
read(1,*) nray(1),nray(2)
close(1)


do iiips=1,2
    if(nray(iiips).eq.0) cycle
    write(ppss,'(i1)')iiips

    write(*,*)' *****************************************'
    if(iiips.eq.1) write(*,*)' P-model:'
    if(iiips.eq.2) write(*,*)' S-model:'

    links=0

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/kod_xy'//ppss//gr//'.dat')
    read(1,*)nxpl,nypl
    allocate(kod_xy(nxpl,nypl,2))
    do iy=1,nypl
        read(1,*)((kod_xy(ix,iy,i2),i2=1,2),ix=1,nxpl)
    end do
    close(1)

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr'//ppss//gr//'.dat')
    nmax=0
    do n=1,nypl
        read(1,*)nt,y
        !write(*,*)nt,y
        if(nt.gt.nmax) nmax=nt
        if(nt.eq.0) cycle
        do i=1,nt
	        read(1,*)x,z,l
        end do
    end do
    close(1)
    write(*,*)' nmax=',nmax

    allocate(xtop(nmax,nypl),ztop(nmax,nypl),popor(nmax,nypl))

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr'//ppss//gr//'.dat')
    do iy=1,nypl
        read(1,*)ntop(iy),ylevel(iy)
        !write(*,*)ntop(iy),ylevel(iy)
        if(ntop(iy).eq.0) cycle
        do inod=1,ntop(iy)
            read(1,*)xtop(inod,iy),ztop(inod,iy),popor(inod,iy)
        end do
    end do
    close(1)

    ilink=0

    ! finding links in X-direction
    do iy=1,nypl
        yy=ylevel(iy)
        if(ntop(iy).eq.0) cycle
        do ix=1,nxpl-1
            if(kod_xy(ix,iy,1).eq.0.or.kod_xy(ix+1,iy,1).eq.0) cycle
            iz11=kod_xy(ix,iy,1)
            iz12=kod_xy(ix,iy,2)
            iz21=kod_xy(ix+1,iy,1)
            iz22=kod_xy(ix+1,iy,2)
            !write(*,*)' iz11, iz12=',iz11, iz12
            !write(*,*)' iz21, iz22=',iz21, iz22
            nz1=(iz12-iz11)+1
            do iz=1,nz1
                izcur=iz-1+iz11
                zver1(iz)=ztop(izcur,iy)
                !write(*,*)' izcur=',izcur,' xtop(izcur,iy)=',xtop(izcur,iy),' ztop(izcur,iy)=',ztop(izcur,iy)
            end do
            nz2=(iz22-iz21)+1
            do iz=1,nz2
                izcur=iz-1+iz21
                zver2(iz)=ztop(izcur,iy)
                    !write(*,*)' izcur=',izcur,' xtop(izcur,iy)=',xtop(izcur,iy),' ztop(izcur,iy)=',ztop(izcur,iy)
            end do
    !        write(*,'(5f7.2)')(zver1(iz),iz=1,nz1)
    !        write(*,'(5f7.2)')(zver2(iz),iz=1,nz2)

            icur1=1
            icur2=1

    81      k_lay_1=icur1-1+iz11
            k_tot_1=popor(k_lay_1,iy)

            k_lay_2=icur2-1+iz21
            k_tot_2=popor(k_lay_2,iy)

            ilink=ilink+1
            links(1,ilink)=k_tot_1
            links(2,ilink)=k_tot_2
            !write(*,*)' ilink=',ilink,' k_tot_1=',k_tot_1,' k_tot_2=',k_tot_2

            if(icur1.lt.nz1) then
                inext1=icur1+1
            else
                inext1=0
            end if

            if(icur2.lt.nz2) then
                inext2=icur2+1
            else
                inext2=0
            end if

            if(inext1*inext2.ne.0) then
                zcur1=zver1(icur1)
                zcur2=zver2(icur2)
                znext1=zver1(inext1)
                znext2=zver2(inext2)
                if(abs(zcur1-znext2).gt.abs(zcur2-znext1)) then
                    icur1=inext1
                else
                    icur2=inext2
                end if
            else if(inext1.ne.0) then
                icur1=inext1
            else if(inext2.ne.0) then
                icur2=inext2
            else
                !pause
                cycle
            end if
            goto 81

        end do
        !write(*,*)' iy=',iy,' yy=',yy,' ilink=',ilink
        !pause
    end do

    write(*,*)' links in X-direction: ilink=',ilink


    ! finding links in Y-direction
    do ix=1,nxpl
        do iy=1,nypl-1
            yy=ylevel(iy)
            if(kod_xy(ix,iy,1).eq.0.or.kod_xy(ix,iy+1,1).eq.0) cycle
            iz11=kod_xy(ix,iy,1)
            iz12=kod_xy(ix,iy,2)
            iz21=kod_xy(ix,iy+1,1)
            iz22=kod_xy(ix,iy+1,2)
            !write(*,*)' iz11, iz12=',iz11, iz12
            !write(*,*)' iz21, iz22=',iz21, iz22
            nz1=(iz12-iz11)+1
            do iz=1,nz1
                izcur=iz-1+iz11
                zver1(iz)=ztop(izcur,iy)
                !write(*,*)' izcur=',izcur,' xtop(izcur,iy)=',xtop(izcur,iy),' ztop(izcur,iy)=',ztop(izcur,iy)
            end do
            nz2=(iz22-iz21)+1
            do iz=1,nz2
                izcur=iz-1+iz21
                zver2(iz)=ztop(izcur,iy+1)
                    !write(*,*)' izcur=',izcur,' xtop(izcur,iy)=',xtop(izcur,iy),' ztop(izcur,iy)=',ztop(izcur,iy)
            end do
            !write(*,'(5f7.2)')(zver1(iz),iz=1,nz1)
            !write(*,'(5f7.2)')(zver2(iz),iz=1,nz2)

            icur1=1
            icur2=1

    82      k_lay_1=icur1-1+iz11
            k_tot_1=popor(k_lay_1,iy)

            k_lay_2=icur2-1+iz21
            k_tot_2=popor(k_lay_2,iy+1)

            ilink=ilink+1
            links(1,ilink)=k_tot_1
            links(2,ilink)=k_tot_2
            !write(*,*)' ilink=',ilink,' k_tot_1=',k_tot_1,' k_tot_2=',k_tot_2

            if(icur1.lt.nz1) then
                inext1=icur1+1
            else
                inext1=0
            end if

            if(icur2.lt.nz2) then
                inext2=icur2+1
            else
                inext2=0
            end if

            if(inext1*inext2.ne.0) then
                zcur1=zver1(icur1)
                zcur2=zver2(icur2)
                znext1=zver1(inext1)
                znext2=zver2(inext2)
                if(abs(zcur1-znext2).gt.abs(zcur2-znext1)) then
                    icur1=inext1
                else
                    icur2=inext2
                end if
            else if(inext1.ne.0) then
                icur1=inext1
            else if(inext2.ne.0) then
                icur2=inext2
            else
                cycle
            end if
            goto 82

        end do
        !write(*,*)' ix=',ix,' ilink=',ilink
    end do

    write(*,*)' links in Y-direction: ilink=',ilink

    open(12,file='../../../DATA/'//ar//'/'//md//'/data/link_hor'//ppss//gr//'.dat',form='binary')
    write(12)ilink
    do i=1,ilink
	    write(12)links(1,i),links(2,i)
    end do
    close(12)

    ! finding links in Z-direction
    ilink=0
    do iy=1,nypl
        yy=ylevel(iy)
        if(ntop(iy).eq.0) cycle
        do ix=1,nxpl-1
            if(kod_xy(ix,iy,1).eq.0) cycle
            iz1=kod_xy(ix,iy,1)
            iz2=kod_xy(ix,iy,2)
            if(iz1.eq.iz2)cycle
            !write(*,*)' ix=',ix,' iy=',iy,' iz1=',iz1,' iz2=',iz2
            do iz=iz1,iz2-1
                k_tot_1=popor(iz,iy)
                k_tot_2=popor(iz+1,iy)

                ilink=ilink+1
                links(1,ilink)=k_tot_1
                links(2,ilink)=k_tot_2
                !write(*,*)' ilink=',ilink,' k_tot_1=',k_tot_1,' k_tot_2=',k_tot_2
            
            end do
        end do
    end do
    write(*,*)' links in Z-direction: ilink=',ilink

    open(12,file='../../../DATA/'//ar//'/'//md//'/data/link_ver'//ppss//gr//'.dat',form='binary')
    write(12)ilink
    do i=1,ilink
	    write(12)links(1,i),links(2,i)
    end do
    close(12)
    deallocate(xtop,ztop,popor,kod_xy)

end do

stop
end