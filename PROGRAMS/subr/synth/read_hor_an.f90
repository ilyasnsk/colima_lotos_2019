subroutine read_hor_an(ar,md)
character*8 ar,md
character*5 fig
real xcur(1000),ycur(1000)

common/center/fi0,tet0
common/pi/pi,per
common/dv_hor_an/nan,anom(40,400,2),val(40,4),zan1(40),zan2(40),nnod(40)
common/keys/key_ft1_xy2

iscale=0
open(2,file='../../../DATA/'//ar//'/'//md//'/forms/scaling.dat')
read(2,*,end=392)fmap1,fmap2
read(2,*)tmap1,tmap2
read(2,*)rotate
iscale=1
392 close(2)

sina=sin(rotate*per)
cosa=cos(rotate*per)


open(1,file='../../../DATA/'//ar//'/'//md//'/anomaly.dat')
read(1,*) 
read(1,*) 
read(1,*) nform



do iform=1,nform
	read(1,*) 
	read(1,*) fig
	read(1,*) xpmin,xpmax,ypmin,ypmax
	read(1,*) val(iform,1),val(iform,2),val(iform,3),val(iform,4)
	read(1,*) zan1(iform),zan2(iform)

	open(2,file='../../../DATA/'//ar//'/'//md//'/forms/'//fig//'.bln')
	read(2,*) ncurve
	xmin=10000
	xmax=-10000
	ymin=10000
	ymax=-10000
	do i=1,ncurve
		read(2,*)xcur(i),ycur(i)
		if(xcur(i).lt.xmin) xmin=xcur(i)
		if(xcur(i).gt.xmax) xmax=xcur(i)
		if(ycur(i).lt.ymin) ymin=ycur(i)
		if(ycur(i).gt.ymax) ymax=ycur(i)
	end do
	close(2)
	!write(*,*)' xmin=',xmin,' xmax=',xmax
	!write(*,*)' zmin=',zmin,' zmax=',zmax
	!pause

	dfdx=(xpmax-xpmin)/(xmax-xmin)
	dtdy=(ypmax-ypmin)/(ymax-ymin)
	if(abs(xpmin-xpmax).lt.0.000001) then
		do i=1,ncurve
			anom(iform,i,1)=xcur(i)
			anom(iform,i,2)=ycur(i)
			!write(*,*)anom(iform,i,1),anom(iform,i,2)
		end do
		nnod(iform)=ncurve+1
		anom(iform,nnod(iform),1)=xcur(1)
		anom(iform,nnod(iform),2)=ycur(1)
	else
		do i=1,ncurve
			anom(iform,i,1)=xpmin+dfdx*(xcur(i)-xmin)
			anom(iform,i,2)=ypmin+dtdy*(ycur(i)-ymin)
		end do
		nnod(iform)=ncurve+1
		anom(iform,nnod(iform),1)=xpmin+dfdx*(xcur(1)-xmin)
		anom(iform,nnod(iform),2)=ypmin+dtdy*(ycur(1)-ymin)
	end if
end do

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
		!write(*,*)' nnod=',nnod(iform)
		do i=1,nnod(iform)
			x1=anom(iform,i,1)
			y1=anom(iform,i,2)

			x2=fmap1+((fmap2-fmap1)/(xmax-xmin))*(x1-xmin)
			y2=tmap1+((tmap2-tmap1)/(ymax-ymin))*(y1-ymin)

			anom(iform,i,1)=x2
			anom(iform,i,2)=y2
		end do
	end do
end if

!write(*,*)' key_ft1_xy2=',key_ft1_xy2
!pause
do iform=1,nform
	!write(*,*)' nnod=',nnod(iform)
	do i=1,nnod(iform)
		fi=anom(iform,i,1)
		tet=anom(iform,i,2)
                if(key_ft1_xy2.eq.1) then
		    call SFDEC(FI,TET,0.,X,Y,Z,fi0,tet0)
                else
                    x=fi
                    y=tet
                end if
		anom(iform,i,1)=x
		anom(iform,i,2)=y
		!write(*,*)x,y
	end do
end do
close(1)




nan=nform
!write(*,*)' nan=',nan

return
end