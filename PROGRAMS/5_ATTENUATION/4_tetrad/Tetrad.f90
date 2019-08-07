real xtop1(9900),ztop1(9900),xtmp(9900),ztmp(9900)
INTEGER otr1(2,9200),otr(2,9200,200),notr(200),ntop(200),nver(200)
integer ssame1(9500),ndoub(200)

real xtop2(9900),ztop2(9900),xtop(9900,200),ztop(9900,200)
real verlin1(990),verlin2(990),verlin(990,200)
real zver1(345,990),zver2(345,990),ylevel(200)
INTEGER otr2(2,9200)
integer tet(4,9000),obsos(310)
integer obr1(335,395), obr2(335,395), nn1(395), nn2(395)
integer sosedi(330,9000)
integer npar(200),nray(2)

real ztest(4),vertet(4,3)
integer iuz(4),nurr(4)

character*3 lv
character*8 ar,md,line
character*1 it,ppss,gr,ps

ioddeven=0
open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)ips
read(1,*)igr
close(1)
write(ps,'(i1)')ips
write(gr,'(i1)')igr


write(*,*)' execution of Tetrad'
write(*,*)' ar=',ar,' md=',md,' ps=',ps,' gr=',gr


nsosmax=100

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

nxpl=(xlim2-xlim1)/dxpl
nypl=(ylim2-ylim1)/dypl
nzpl=(zlim2-zlim1)/dzpl




open(1,file='../../../DATA/'//ar//'/'//md//'/data/lev_att'//gr//ps//'.dat')
i=0

722 i=i+1
read(1,*,end=721)npar(i),ylevel(i)
!write(*,*)npar(i),ylevel
goto 722
721 nypl=i-1
close(1)

!write(*,*)' nypl=',nypl


open(3,file='../../../TMP_files/tmp/ver_att'//gr//ps//'.dat')
open(11,file='../../../TMP_files/tmp/otr_att'//gr//ps//'.dat')

open(2,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
do nur=1,nypl
read(2,*)ntop(nur)
!	write(*,*)ntop(nur)
do i=1,ntop(nur)
	read(2,*)xtop(i,nur),ztop(i,nur),ipop
!		write(*,*)xtop(i,nur),ztop(i,nur),ipop
end do
		
do i=2,ntop(nur)
	if(abs(ztop(i-1,nur)-ztop(i,nur)).lt.0.1) then
		ztop(i-1,nur)=ztop(i-1,nur)-0.5
		ztop(i,nur)=ztop(i,nur)+0.5
	end if
end do
!write(*,*)' nur=',nur,' notr=',notr(nur),' ntop=',ntop(nur),' nver=',nver(nur)
end do
close(2)




! Compute otrezki at each level

do n=1,nypl
!do n=32,32

    ntop1=ntop(n)
    do i=1,ntop1
        xtop1(i)=xtop(i,n)
        ztop1(i)=ztop(i,n)
    end do
    !	write(*,*)' ntop1=',ntop1
    !	do i=1,ntop1
    !		write(*,*)xtop1(i),ztop1(i)
    !	end do
    otr1=0
    notr1=0
    npoinonver=1
    iverline=1
    zver1(1,1)=ztop1(1)
    obr1(1,1)=1
    do l = 2, ntop1
        if (abs(xtop1(l)-xtop1(l-1)).gt.0.0001) then
        !	write(*,*)'    iverline=',iverline
	        nn1(iverline)=npoinonver
	        iverline=iverline+1
	        npoinonver=0
        endif
        npoinonver=npoinonver+1
        zver1(npoinonver,iverline)=ztop1(l)
        obr1(npoinonver,iverline)=l
    enddo
    nn1(iverline)=npoinonver

    nvert1=iverline


    do k=1,nvert1
        do j=1,(nn1(k)-1)                 
            notr1=notr1+1
            otr1(1,notr1)=obr1(j,k)
            otr1(2,notr1)=obr1(j+1,k)
        end do
    end do
    !	write(*,*)' notr1=',notr1

    do iver=1,nvert1-1
        i1=1
        i2=1
        notr1=notr1+1
        otr1(1,notr1)=obr1(1,iver)
        otr1(2,notr1)=obr1(1,iver+1)
        !		write(*,'(10f7.1)')(zver1(i,iver),i=1,nn1(iver))
        !		write(*,'(10f7.1)')(zver1(i,iver+1),i=1,nn1(iver+1))
        1		continue
        !		write(*,*)' i1=',i1,' i2=',i2
        !		pause
        if(i1.ge.nn1(iver).and.i2.ge.nn1(iver+1))cycle
        if(i1.ge.nn1(iver))goto 332
        if(i2.ge.nn1(iver+1))goto 333

	! Double points:

        id1=1
        id2=1
        id11=1
        id22=1
        id01=1
        id02=1
        if(i1+2.le.nn1(iver)) then
	        if (abs(zver1(i1,iver)-zver1(i1+1,iver)).lt.1.11) id1=2
	        if (abs(zver1(i1+2,iver)-zver1(i1+1,iver)).lt.1.11) id11=2
        !			write(*,*)' id1=',id1,' id11=',id11
        end if
        if(i2+2.le.nn1(iver+1)) then
	        if (abs(zver1(i2,iver+1)-zver1(i2+1,iver+1)).lt.1.11) id2=2
	        if (abs(zver1(i2+2,iver+1)-zver1(i2+1,iver+1)).lt.1.11) id22=2
        !			write(*,*)' id2=',id2,' id22=',id22
        end if
        if(i1.ne.1) then
	        if (abs(zver1(i1,iver)-zver1(i1-1,iver)).lt.1.11) id01=2
        end if
        if(i2.ne.1) then
	        if (abs(zver1(i2,iver+1)-zver1(i2-1,iver+1)).lt.1.11) id02=2
        end if
        331		continue

        if(id1.eq.2.and.id2.eq.2) then
	        notr1=notr1+1
	        otr1(1,notr1)=obr1(i2+1,iver+1)
	        otr1(2,notr1)=obr1(i1+1,iver)
	        i1=i1+1
	        i2=i2+1
	        goto 1
        end if

        if(id1.eq.2.and.id2.eq.1) goto 332
        if(id2.eq.2.and.id1.eq.1) goto 333

        dz1=abs(zver1(i1,iver)-zver1(i2+1,iver+1))
        dz2=abs(zver1(i2,iver+1)-zver1(i1+1,iver))
        if(dz1*10.lt.dz2) goto 332
        if(dz2*10.lt.dz1) goto 333



        if (dz1.lt.dz2) goto 332

        !*******************************************************************
        333		notr1=notr1+1
        otr1(1,notr1)=obr1(i2,iver+1)
        otr1(2,notr1)=obr1(i1+1,iver)
        i1=i1+1
        goto 1
        !*******************************************************************
        332		notr1=notr1+1
        otr1(1,notr1)=obr1(i1,iver)
        otr1(2,notr1)=obr1(i2+1,iver+1)
        i2=i2+1
        goto 1
        !*******************************************************************
        !		write(*,*)' notr1=',notr1
    end do
    !write(*,*)n,' notr1=',notr1
    write(11,*) notr1
    write(11,*) ((otr1(i,j), i=1,2), j=1,notr1)

    yy=ylevel(n)
    !call decsf(0.,yy,0.,fi0,tet0,fff,tetsec,h)
    !itetsec=int(tetsec*100.)
    !write(lv,'(i3)')n
    !open(12,file='../../../TMP_files/otrezki/otr'//lv//ppss//gr//'.bln')
    !do i=1,notr1
	    !x11=xtop1(otr1(1,i))
	    !z11=ztop1(otr1(1,i))
	    !x22=xtop1(otr1(2,i))
	    !z22=ztop1(otr1(2,i))
	    !if(n.eq.32)write(*,'(3i5,4f8.3)')i,otr1(1,i),otr1(2,i),x11,x22,z11,z22
	    !write(*,*)i,x11,x22,z11,z22
	    !write(12,*)2
	    !write(12,*)x11,-z11
	    !write(12,*)x22,-z22
    !end do
    !close(12)

    !	do j=1,notr1
    !		write(*,*) (otr1(i,j), i=1,2)
    !	end do

    ssame1=0
    isos=0.
    do i=1,ntop1
        isos=isos+1
        ssame1(isos)=0
        isos=isos+1
        ssame1(isos)=i
        do  l= 1, notr1
            if (otr1(1,l).eq.i) then
                isos=isos+1
                ssame1(isos)=otr1(2,l)
            endif
            if (otr1(2,l).eq.i) then
                isos=isos+1
                ssame1(isos)=otr1(1,l)
            endif
        enddo
    end do
end do
close(11)
close(2)
close(3)
!pause

open(3,file='../../../TMP_files/tmp/ver_att'//gr//ps//'.dat')
open(4,file='../../../TMP_files/tmp/otr_att'//gr//ps//'.dat')
open(12,file='../../../TMP_files/tmp/tetr_att'//gr//ps//'.dat')

do nur=1,nypl
    read(4,*) notr(nur)
    read(4,*) ((otr(i,j,nur), i=1,2), j=1,notr(nur))
    read(3,*)nver(nur)
    read(3,*)(verlin(i,nur),i=1,nver(nur))
end do

do nur=1,nypl-1

    nur1=nur
    nur2=nur+1

    notr1=notr(nur1) 
    notr2=notr(nur2) 



    !	write(*,*)' nur1=',nur1,' nur2=',nur2
    !	write(*,*)' notr1=',notr1,' notr2=',notr2

    do j=1,2
	    do i=1,notr1
		    otr1(j,i)=otr(j,i,nur1)
	    end do
	    do i=1,notr2
		    otr2(j,i)=otr(j,i,nur2)
	    end do
    end do
    ntop1=ntop(nur1) 
    ntop2=ntop(nur2) 
    !	write(*,*)' ntop1=',ntop1,' ntop2=',ntop2
    do i=1,ntop1
        xtop1(i)=xtop(i,nur1)
        ztop1(i)=ztop(i,nur1)
    end do
    do i=1,ntop2
        xtop2(i)=xtop(i,nur2)
        ztop2(i)=ztop(i,nur2)
    end do
    nver1=nver(nur1) 
    nver2=nver(nur2) 
    !	write(*,*)' nver1=',nver1,' nver2=',nver2
    !	pause
    ! compute SOSEDI at first level


    npoinonver=1
    iverline=1
    zver1(1,1)=ztop1(1)
    obr1(1,1)=1
    do l = 2, ntop1
        if (abs(xtop1(l)-xtop1(l-1)).gt.0.0001) then
            !	write(*,*)'    iverline=',iverline
            nn1(iverline)=npoinonver
            !write(*,'(20i4)')(obr1(i,iverline),i=1,nn1(iverline))
            iverline=iverline+1
            npoinonver=0
        endif
        npoinonver=npoinonver+1
        zver1(npoinonver,iverline)=ztop1(l)
        obr1(npoinonver,iverline)=l
    enddo
    nn1(iverline)=npoinonver
    nvert1=iverline
    !	write(*,*)

    npoinonver=1
    iverline=1
    zver2(1,1)=ztop2(1)
    obr2(1,1)=1
    do l = 2, ntop2
        if (abs(xtop2(l)-xtop2(l-1)).gt.0.0001) then
    !	write(*,*)'    iverline=',iverline
	    nn2(iverline)=npoinonver
    !write(*,'(20i4)')(obr2(i,iverline),i=1,nn2(iverline))
	    iverline=iverline+1
	    npoinonver=0
        endif
        npoinonver=npoinonver+1
        zver2(npoinonver,iverline)=ztop2(l)
        obr2(npoinonver,iverline)=l
    enddo
    nn2(iverline)=npoinonver
    nvert2=iverline

    sosedi=0
    do i=1,ntop1
        nsos=0
        do k = 1, notr1
            if (otr1(1,k).eq.i) then        
	            nsos=nsos+1                  
	            sosedi(nsos,i)=otr1(2,k)     
            endif
            if (otr1(2,k).eq.i) then
	            nsos=nsos+1
	            sosedi(nsos,i)=otr1(1,k)
            endif
        enddo
        !		write(*,'(20i4)')i,(sosedi(k,i),k=1,nsos)
    enddo
    !	do k=1,notr1
    !		write(*,*)otr1(1,k),otr1(2,k)
    !	end do
    !	write(*,'(20i4)')111,(sosedi(k,111),k=1,10)


    ! compute SOSEDI at second level
    do i=1,ntop2
        nsos=0
        do k = 1, notr2
            if (otr2(1,k).eq.i) then
	            nsos=nsos+1
	            sosedi(nsos,i+ntop1)=otr2(2,k)+ntop1
            endif
            if (otr2(2,k).eq.i) then
	            nsos=nsos+1
	            sosedi(nsos,i+ntop1)=otr2(1,k)+ntop1
            endif
        enddo
        !		write(*,'(20i4)')i,(sosedi(k,i+ntop1),k=1,nsos)
    enddo



    tet=0
    ntet=0
    npver=1                             
    zver1(1,1)=ztop1(1)
    obr1(1,1)=1
    k=1                                 
    verlin1(k)=xtop1(1)
    do m=2,ntop1                        
        if (abs(xtop1(m)-xtop1(m-1)).gt.0.0001)then
            nn1(k)=npver
            !			write(*,*)' k=',k,' npver111=',npver
            k=k+1
            verlin1(k)=xtop1(m)
            npver=0
        end if
        npver=npver+1
        zver1(npver,k)=ztop1(m)
        obr1(npver,k)=m
    end do
    nn1(k)=npver
    nver1=k

    npver=1                            
    zver2(1,1)=ztop2(1)
    obr2(1,1)=1+ntop1
    k=1
    verlin2(k)=xtop2(1)
    do m=2,ntop2
        if (abs(xtop2(m)-xtop2(m-1)).gt.0.0001)then
            nn2(k)=npver
            !	write(*,'(20i4)')(obr2(i,k)-ntop1,i=1,nn2(k))
            !			write(*,*)' k=',k,' npver222=',npver
            k=k+1
            verlin2(k)=xtop2(m)
            npver=0
        end if
        npver=npver+1
        zver2(npver,k)=ztop2(m)
        obr2(npver,k)=m+ntop1
    end do
    nn2(k)=npver
    nver2=k

    ivert1=1
    ivert2=1
    markver=1

    ivertot=0

    10  x1=verlin1(ivert1)
    x2=verlin2(ivert2)
    if(markver.eq.1) goto 229
    if (ivert1.eq.nver1) then
	    x11=1.e10
    else
	    x11=verlin1(ivert1+1)
    endif

    if (ivert2.eq.nver2) then
	    x22=1.e10
    else
	    x22=verlin2(ivert2+1)
    endif

    if(ivert1.eq.nver1.and.ivert2.eq.nver2) goto 15

    if (abs(x11-x2).gt.abs(x1-x22)) then
	    ivert2=ivert2+1
	    nurtet=1
    else
	    ivert1=ivert1+1
	    nurtet=0
    endif    



    229 markver=0
    !	write(*,*)' ivert1=',ivert1,' ivert2=',ivert2
    !	write(*,*)' verli1=',verlin1(ivert1),' verli2=',verlin2(ivert2)
    ivertot=ivertot+1
    yris1=-(ivertot-1)*10-2
    yris2=-ivertot*10+2

    !	do i=1,nn1(ivert1)
    !		write(32,*)zver1(i,ivert1),yris1
    !	end do
    !	do i=1,nn2(ivert2)
    !		write(32,*)zver2(i,ivert2),yris2
    !	end do


    i1=obr1(1,ivert1)
    i2=obr2(1,ivert2)

    !	write(*,'(20i4)')(obr1(i,ivert1),i=1,nn1(ivert1))
    !	write(*,'(20i4)')(obr2(i,ivert2)-ntop1,i=1,nn2(ivert2))
    !	write(*,'(10f7.1)')(zver1(i,ivert1),i=1,nn1(ivert1))
    !	write(*,'(10f7.1)')(zver2(i,ivert2),i=1,nn2(ivert2))
    !	write(*,*)
    !	write(*,*)' i1=',i1,' i2=',i2-ntop1
    !	write(*,*)ztop1(i1),ztop2(i2-ntop1)
    !	write(31,*)2
    !	write(31,*)ztop1(i1),yris1
    !	write(31,*)ztop2(i2-ntop1),yris2
    !	pause

    include 'otz.fl'   
    include 'tetr.fl'
    !	write(*,*)' ntet=',ntet
    !	if(ntet.ne.0)write(*,*)ntet,(tet(i,ntet),i=1,4)

    l1=1
    l2=1

    13  l11=l1+1
    l22=l2+1

    id01=1
    id1=1
    id11=1
    if(l1+2.le.nn1(ivert1)) then
	    if (abs(zver1(l1,ivert1)-zver1(l1+1,ivert1)).lt.1.11) id1=2
	    if (abs(zver1(l1+1,ivert1)-zver1(l1+2,ivert1)).lt.1.11) id11=2
    !			write(*,*)' id1=',id1,' id11=',id11
    end if
    if(l1.ne.1) then
	    if (abs(zver1(l1,ivert1)-zver1(l1-1,ivert1)).lt.1.11) id01=2
    end if

    id02=1
    id2=1
    id22=1
    if(l2+2.le.nn2(ivert2)) then
	    if (abs(zver2(l2,ivert2)-zver2(l2+1,ivert2)).lt.1.11) id2=2
	    if (abs(zver2(l2+1,ivert2)-zver2(l2+2,ivert2)).lt.1.11) id22=2
    !			write(*,*)' id1=',id1,' id11=',id11
    end if
    if(l2.ne.1) then
	    if (abs(zver1(l2,ivert2)-zver1(l2-1,ivert1)).lt.1.11) id02=2
    end if

    !write(*,*)' id01=',id01,' id1=',id1,' id11=',id11
    !write(*,*)' id02=',id02,' id2=',id2,' id22=',id22
    if(l1.eq.nn1(ivert1).and.l2.eq.nn2(ivert2)) goto 14

    if(l1.eq.nn1(ivert1))goto 772
    if(l2.eq.nn2(ivert2))goto 771

    if(id1.eq.2.and.id2.eq.2) then
	    l1=l11
	    l2=l22
	    goto 773
    end if

    if(id1.eq.2.and.id2.eq.1) goto 772
    if(id2.eq.2.and.id1.eq.1) goto 771
    dr1=abs(zver2(l2,ivert2)-zver1(l11,ivert1))
    dr2=abs(zver2(l22,ivert2)-zver1(l1,ivert1))

    if(dr1*10.lt.dr2) goto 771
    if(dr2*10.lt.dr1) goto 772





    !	write(*,*)' z1=',zver1(l1,ivert1),' z22=',zver2(l22,ivert2)
    !	write(*,*)' dr2=',dr2,' l1=',l1,' l22=',l22
    !	write(*,*)' verlin1=',verlin1(ivert1),' verlin2=',verlin2(ivert2)
    !	pause


    if (dr2.lt.dr1) goto 772

    771	continue
    l1=l11
    goto 773
    772	l2=l22
    773 continue
    !	write(*,*)' l1=',l1,' l2=',l2

    i1=obr1(l1,ivert1)
    i2=obr2(l2,ivert2)
    !	write(31,*)2
    !	write(31,*)ztop1(i1),yris1
    !	write(31,*)ztop2(i2-ntop1),yris2
    !	write(*,*)' i1=',i1,' i2=',i2-ntop1
    x1=xtop1(i1)
    z1=ztop1(i1)
    !	write(*,*)' x1=',x1,' z1=',z1
    x2=xtop2(i2-ntop1)
    z2=ztop2(i2-ntop1)
    !	write(*,*)' x2=',x2,' z2=',z2
    !	pause
    !	write(22,*)2
    !	write(22,*)x1,-z1
    !	write(22,*)x2,-z2
    !	pause
    do isos=1,nsosmax
	    if(sosedi(isos,i1).eq.0) exit
    end do
    sosedi(isos,i1)=i2
    !	write(*,*)' isos=',isos
    do isos=1,nsosmax
	    if(sosedi(isos,i2).eq.0) exit
    end do
    sosedi(isos,i2)=i1
    !	write(*,*)' isos=',isos
    nold=ntet

    include'tetr.fl'
    !		write(*,'(i5,5x,4i5)')ntet,(tet(i,ntet),i=1,4)

    !	if(nold.ne.ntet)then
    !		if(ivert1.eq.11.and.ivert2.eq.17) then
    !			do nn=nold+1,ntet
    !			end do
    !		end if
    !	end if
    !	write(*,*)' ntet=',ntet
    !	if(ntet.ne.0)write(*,*)ntet,(tet(i,ntet),i=1,4)
    goto 13

    14  continue



    goto 10

    15  continue
    write(12,*) ntet
    if(mod(nur,20).eq.0)write(*,*)' nur=',nur,' ntet=',ntet,' ntop=',ntop(nur1)
    do n=1,ntet
        !		write(*,*)(tet(i,n),i=1,4)
        write(12,*)(tet(i,n),i=1,4)
        do i1=1,4
            itt1=tet(i1,n)
            do i2=1,4
                if(i2.eq.i1) cycle
                itt2=tet(i2,n)
                if(itt1.le.ntop1.and.itt2.gt.ntop1)then
                    !					write(21,*)2
                    x1=xtop1(itt1)
                    z1=ztop1(itt1)
                    !					write(21,*)x1,-z1
                    x2=xtop2(itt2-ntop1)
                    z2=ztop2(itt2-ntop1)
                    !					write(21,*)x2,-z2
                else if(itt1.le.ntop1.and.itt2.le.ntop1)then
                    !					write(22,*)2
                    x1=xtop1(itt1)
                    z1=ztop1(itt1)
                    !					write(22,*)x1,-z1
                    x2=xtop1(itt2)
                    !					z2=ztop1(itt2)
                    !					write(22,*)x2,-z2
                else if(itt1.gt.ntop1.and.itt2.gt.ntop1)then
                    !					write(23,*)2
                    x1=xtop2(itt1-ntop1)
                    z1=ztop2(itt1-ntop1)
                    !					write(23,*)x1,-z1
                    x2=xtop2(itt2-ntop1)
                    z2=ztop2(itt2-ntop1)
                    !					write(23,*)x2,-z2
                end if
            end do
        end do
    end do
987    continue
end do


stop
end
