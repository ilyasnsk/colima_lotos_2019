subroutine read_ini_model(ar,md)
character*8 ar,md
character*5 fig
real xcur(1000),zcur(1000)
common/center/fi0,tet0



common/dv_an/nan,anom(40,400,2),val_p(40),val_s(40),yan1(40),yan2(40),nnod(40)
common/profile/xa0(40),ya0(40),xb0(40),yb0(40)

nform=0
open(1,file='../../../DATA/'//ar//'/'//md//'/ini_model.dat',status='old',err=331)
read(1,*,end=331) 
read(1,*) 
read(1,*) nform



do iform=1,nform
	read(1,*) 
	read(1,*) fia,teta,fib,tetb
	read(1,*) fig
	read(1,*) xpmin,xpmax,zpmin,zpmax
	read(1,*) val_p(iform),val_s(iform)
	read(1,*) yan1(iform),yan2(iform)

	call SFDEC(fia,teta,0.,xa,ya,Z,fi0,tet0)
	call SFDEC(fib,tetb,0.,xb,yb,Z,fi0,tet0)

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

331 nan=nform
!write(*,*)' nan=',nan

return
end