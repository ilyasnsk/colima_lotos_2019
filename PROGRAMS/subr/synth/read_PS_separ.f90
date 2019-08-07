subroutine read_PS_separ(ar,md)
character*8 ar,md
character*5 fig
real xcur(1000),zcur(1000)
common/center/fi0,tet0
common/pi/pi,per



common/dv_PS_separ/nan,anom(40,400,2),ips_separ(40),val_separ(40),yan1(40),yan2(40),nnod(40)
common/profile/xa0(40),ya0(40),xb0(40),yb0(40)


iscale=0
open(2,file='../../../DATA/'//ar//'/'//md//'/forms/scaling.dat')
read(2,*,end=392)xmap1,xmap2
read(2,*)ymap1,ymap2
read(2,*)rotate
iscale=1
392 close(2)

sina=sin(rotate*per)
cosa=cos(rotate*per)




nform=0
open(1,file='../../../DATA/'//ar//'/'//md//'/anomaly.dat')
read(1,*,end=331) 
read(1,*) 
read(1,*) nform



do iform=1,nform
	read(1,*) 
	read(1,*) ips_separ(iform)
	read(1,*) fia,teta,fib,tetb
	read(1,*) fig
	read(1,*) xpmin,xpmax,zpmin,zpmax
	read(1,*) val_separ(iform)
	read(1,*) yan1(iform),yan2(iform)

	if(fi0.lt.-900) then
		xa=fia
		ya=teta
		xb=fib
		yb=tetb
	else
		call SFDEC(fia,teta,0.,xa,ya,Z,fi0,tet0)
		call SFDEC(fib,tetb,0.,xb,yb,Z,fi0,tet0)
	end if


	!write(*,*)' xa=',xa,' ya=',ya
	!write(*,*)' xb=',xb,' yb=',yb
	xa0(iform)=xa
	xb0(iform)=xb
	ya0(iform)=ya
	yb0(iform)=yb

	open(2,file='../../../DATA/'//ar//'/'//md//'/forms/'//fig//'.bln')
	read(2,*) ncurve
	xmin=10000
	xmax=-10000
	zmin=10000
	zmax=-10000
	do i=1,ncurve
		read(2,*)xcur(i),zcur(i)
		!write(*,*) xcur(i),zcur(i)
		if(xcur(i).lt.xmin) xmin=xcur(i)
		if(xcur(i).gt.xmax) xmax=xcur(i)
		if(zcur(i).lt.zmin) zmin=zcur(i)
		if(zcur(i).gt.zmax) zmax=zcur(i)
	end do
	close(2)
	!write(*,*)' xmin=',xmin,' xmax=',xmax
	!write(*,*)' zmin=',zmin,' zmax=',zmax
	!pause

	dfdx=(xpmax-xpmin)/(xmax-xmin)
	dtdy=(zpmax-zpmin)/(zmax-zmin)
	if(abs(xpmin-xpmax).lt.0.000001) then
		do i=1,ncurve
			anom(iform,i,1)=xcur(i)
			anom(iform,i,2)=-zcur(i)
			!write(*,*)anom(iform,i,1),anom(iform,i,2)
		end do
		nnod(iform)=ncurve+1
		anom(iform,nnod(iform),1)=xcur(1)
		anom(iform,nnod(iform),2)=-zcur(1)
		!pause
	else
		do i=1,ncurve
			anom(iform,i,1)=xpmin+dfdx*(xcur(i)-xmin)
			anom(iform,i,2)=zpmin+dtdy*(zcur(i)-zmin)
		end do
		nnod(iform)=ncurve+1
		anom(iform,nnod(iform),1)=xpmin+dfdx*(xcur(1)-xmin)
		anom(iform,nnod(iform),2)=zpmin+dtdy*(zcur(1)-zmin)
	end if

	!write(*,*)' nnod=',nnod(iform)
end do
close(1)

if(iscale.eq.1) then
	xmin=99999
	xmax=-99999
	ymin=99999
	ymax=-99999
	do iform=1,nform
		!write(*,*)' nnod=',nnod(iform)
		do i=1,nnod(iform)
			x1=anom(iform,i,1)
			y1=anom(iform,i,2)
			x2=x1*cosa-y1*sina
			y2=x1*sina+y1*cosa
			anom(iform,i,1)=x2
			anom(iform,i,2)=y2
			if(x2.lt.xmin) xmin=x2
			if(x2.gt.xmax) xmax=x2
			if(y2.lt.ymin) ymin=y2
			if(y2.gt.ymax) ymax=y2
		end do
	end do
	!write(*,*)' xmax=',xmax,' xmin=',xmin
	!write(*,*)' ymax=',ymax,' ymin=',ymin

	do iform=1,nform
!		write(*,*)' nnod=',nnod(iform)
		do i=1,nnod(iform)
			x1=anom(iform,i,1)
			y1=anom(iform,i,2)
			!write(*,*)' x1=',x1,' y1=',y1

			x2=xmap1+((xmap2-xmap1)/(xmax-xmin))*(x1-xmin)
			y2=ymap1+((ymap2-ymap1)/(ymax-ymin))*(y1-ymin)
			!write(*,*)' x2=',x2,' y2=',y2

			anom(iform,i,1)=x2
			anom(iform,i,2)=y2
		end do
	end do
end if


331 nan=nform
!write(*,*)' nan=',nan

return
end