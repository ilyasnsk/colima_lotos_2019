character*8 ar,md
character*1 it,gr,ps
integer ntop(200),ntetr(200),tet(:,:,:)
real  xtop(:,:), ztop(:,:),ylevel(200)
integer popor(:,:),notr(200)

integer otr(2,10000,200)

allocatable xtop,ztop,popor,tet


ioddeven=0

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)ips
read(1,*)igr
close(1)

write(*,*)' execution of sosedi'
write(*,*)' ar=',ar,' md=',md,' ps=',ps,' gr=',igr


itdvmod=itstep
!itstep=3
write(ps,'(i1)')ips
write(gr,'(i1)')igr




nlev=nlev+2


open(1,file='../../../DATA/'//ar//'/'//md//'/data/lev_att'//gr//ps//'.dat')
i=0
722 i=i+1
    read(1,*,end=721)ntop(i),ylevel(i)
    !write(*,*)ntop(i),ylevel(i)
    goto 722
721 nypl=i-1
nlev=nypl

!write(*,*)' nypl=',nypl
!pause


open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
nmax=0
do n=1,nypl
    read(1,*)nt
    if(nt.gt.nmax) nmax=nt
    do i=1,ntop(n)
	    read(1,*)x,z,l
    end do
end do
close(1)
!write(*,*)' nmax=',nmax

allocate(xtop(nmax,nypl),ztop(nmax,nypl),popor(nmax,nypl))

open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
do n=1,nypl
    read(1,*)nt
    !	write(*,*)n,' ntop=',nt,ntop(n)
    ntop(n)=nt
    do i=1,ntop(n)
	    read(1,*)xtop(i,n),ztop(i,n),popor(i,n)
    end do
end do
close(1)


open(1,file='../../../TMP_files/tmp/tetr_att'//gr//ps//'.dat')
nmax=0
do n=1,nypl-1
    read(1,*)nt
    if(nt.gt.nmax) nmax=nt
    !write(*,*)n,nt,ntop(n)
    read(1,*)((t,i=1,4),j=1,nt)
end do
close(1)

allocate(tet(4,nmax,nypl-1))

open(1,file='../../../TMP_files/tmp/tetr_att'//gr//ps//'.dat')
do n=1,nypl-1
    if(ntop(n).eq.0) cycle
    read(1,*)ntetr(n)
    read(1,*)((tet(i,j,n),i=1,4),j=1,ntetr(n))
end do
close(1)

notot=0
open(12,file='../../../TMP_files/tmp/sos_att'//gr//ps//'.dat',form='unformatted')
do ilev=1,nlev-1
    notr(ilev)=0
    do itet=1,ntetr(ilev)
        !		write(*,*)itet,(tet(itop1,itet,ilev),itop1=1,4)
        do itop1=1,3
            iuz1=tet(itop1,itet,ilev)
            !			write(*,*)' iuz1=',iuz1,' ntop(ilev)=',ntop(ilev)
            if (iuz1.gt.ntop(ilev)) then
                iuz1=iuz1-ntop(ilev)
                nur1=ilev+1
            else
                nur1=ilev
            end if
            !			write(*,*)' nur1=',nur1,' iuz1=',iuz1
            ipar1=popor(iuz1,nur1)
            if(ipar1.eq.0) cycle
            !			write(*,*)' ipar1=',ipar1
            do itop2=itop1+1,4
                iuz2=tet(itop2,itet,ilev)
                if (iuz2.gt.ntop(ilev)) then
	                iuz2=iuz2-ntop(ilev)
	                nur2=ilev+1
                else
	                nur2=ilev
                end if
                ipar2=popor(iuz2,nur2)
                if(ipar2.eq.0) cycle
                if(notr(ilev).eq.0) goto 777
                do iotr=1,notr(ilev)
	                io1=otr(1,iotr,ilev)
	                io2=otr(2,iotr,ilev)
	                if(io1.eq.ipar1.and.io2.eq.ipar2) goto 778
	                if(io1.eq.ipar2.and.io2.eq.ipar1) goto 778
                end do
                if(ilev.eq.1) goto 777
                do iotr=1,notr(ilev-1)
	                io1=otr(1,iotr,ilev-1)
	                io2=otr(2,iotr,ilev-1)
	                if(io1.eq.ipar1.and.io2.eq.ipar2) goto 778
	                if(io1.eq.ipar2.and.io2.eq.ipar1) goto 778
                end do
                777				continue
                notr(ilev)=notr(ilev)+1
                otr(1,notr(ilev),ilev)=ipar1
                otr(2,notr(ilev),ilev)=ipar2
                x1=xtop(iuz1,nur1)
                z1=ztop(iuz1,nur1)
                y1=ylevel(nur1)
                x2=xtop(iuz2,nur2)
                z2=ztop(iuz2,nur2)
                y2=ylevel(nur2)
                dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
                write(12)ipar1,ipar2,dist
                notot=notot+1
                !if(mod(notot,1000).eq.0)write(*,*)notot,ipar1,ipar2,dist
                778				continue
            end do
        end do
    end do
    if(mod(ilev,20).eq.0) write(*,*)' ilev=',ilev,' notr(ilev)=',notr(ilev),' total=',notot
end do
close(12)

open(12,file='../../../TMP_files/tmp/num_sos_att'//gr//ps//'.dat')
write(12,*)notot
close(12)

deallocate(xtop,ztop,popor,tet)


stop
end