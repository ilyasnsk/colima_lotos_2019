subroutine dv_1_gr_xyz_v(xxx,yyy,zcur,smaxx,   dv,umn)
real wgt(2),dvv(2),xt(3),zt(3),dvt(3)
real amatr(3,3),bmatr(3),sol(3)
integer nonr(540),popr(540,560)


common/pi/pi,per
common/center/fi0,tet0
common/grid/nornt,ornt(4),sinal,cosal,&
		nlev,ylev(300),ntop(300),nobr,&
		xtop(20000,300), ztop(20000,300),n_pop(20000,300),&
		dv_mod(100000),vab_mod(100000)


dv=0
umn=0

!write(*,*)' xxx=',xxx,' yyy=',yyy

xcur=xxx*cosal+yyy*sinal
ycur=-xxx*sinal+yyy*cosal

!write(*,*)' sinal=',sinal,' cosal=',cosal
!write(*,*)' xcur=',xcur,' ycur=',ycur

if(ycur.lt.ylev(1)) return
if(ycur.gt.ylev(nlev)) return

!do ilev=2,nlev
	!write(*,*)ilev,' ntop(ile)=',ntop(ilev)
!end do

do ilev=2,nlev
	y1=ylev(ilev-1)
	y2=ylev(ilev)
	if((ycur-y1)*(ycur-y2).le.0.) exit
end do

!write(*,*)' ilev=',ilev

!write(*,*)' fff=',fff,' ttt=',ttt
!write(*,*)' fi0=',fi0,' tet0=',tet0


!write(*,*)' sinal=',sinal,' cosal=',cosal

!write(*,*)' xcur=',xcur,' ycur=',ycur

xl=xtop(1,ilev)
xr=xtop(ntop(ilev),ilev)
zb=ztop(1,ilev)
zu=ztop(ntop(ilev),ilev)
!write(*,*)xl,xr
!write(*,*)yb,yu
if((xcur-xr)*(xcur-xl).gt.0) return
if((zcur-zb)*(zcur-zu).gt.0) return

dvv=0
wgt=0
do il=1,2
	ile=ilev
	if(il.eq.1) ile=ilev-1

	!write(*,*)' ile=',ile,' ntop(ile)=',ntop(ile)

	wgt(il)=0

	smin=9999999.
	do i=1,ntop(ile)
		if(n_pop(i,ile).eq.0) cycle
		x=xtop(i,ile)
		z=ztop(i,ile)
		y=ylev(ile)
		!write(*,*)' z=',z,' zzz=',zzz
		ss=sqrt((x-xcur)*(x-xcur)+(y-ycur)*(y-ycur)+(z-zcur)*(z-zcur))
		if(ss.lt.smin) smin=ss
		if(ss.lt.smaxx) then
			umnn=1
			goto 71
		end if
	end do
!write(*,*)' smin=',smin,' smaxx=',smaxx


	if(smin.gt.smaxx*2.) then
		umnn=0.
		cycle
	else 
		umnn=(smaxx*2.-smin)/smaxx
	end if
71	continue
!write(*,*)' smin=',smin,' smaxx=',smaxx
!write(*,*)' umn=',umn
	wgt(il)=umnn


	nrow=0
	xold=-100000.
	do itop=1,ntop(ile)
		if(abs(xtop(itop,ile)-xold).gt.0.00001)then
			nrow=nrow+1
			nonr(nrow)=0
			xold=xtop(itop,ile)
		end if
		nonr(nrow)=nonr(nrow)+1
		popr(nrow,nonr(nrow))=itop
	end do
!write(*,*)' nrow=',nrow,' ntop(ile)=',ntop(ile)

	do irow=1,nrow-1
		x1=xtop(popr(irow,1),ile)
		x2=xtop(popr(irow+1,1),ile)
		if((xcur-x1)*(xcur-x2).le.0) exit
	end do
	do iz1=1,nonr(irow)-1
		za1=ztop(popr(irow,iz1),ile)
		zb1=ztop(popr(irow,iz1+1),ile)
		if((zcur-za1)*(zcur-zb1).le.0) goto 897
	end do
	umnn=0.
	cycle

897	continue

	ipop1=n_pop(popr(irow,iz1),ile)
	ipop2=n_pop(popr(irow,iz1+1),ile)
	va1=0
	vb1=0
	if(ipop1.ne.0.and.irow.ne.1)va1=dv_mod(ipop1)
	if(ipop2.ne.0.and.irow.ne.1)vb1=dv_mod(ipop2)

!write(*,*)' xcur=',xcur,' x1=',x1,' x2=',x2
!write(*,*)' zcur=',zcur,' za1=',za1,' zb1=',zb1
!write(*,*)' va1=',va1,' vb1=',vb1

	do iz2=1,nonr(irow+1)-1
		za2=ztop(popr(irow+1,iz2),ile)
		zb2=ztop(popr(irow+1,iz2+1),ile)
		if((zcur-za2)*(zcur-zb2).le.0) goto 898
	end do
	umnn=0.
	cycle

898	continue
	ipop1=n_pop(popr(irow+1,iz2),ile)
	ipop2=n_pop(popr(irow+1,iz2+1),ile)
	va2=0
	vb2=0
	if(ipop1.ne.0.and.irow.ne.nrow-1)va2=dv_mod(ipop1)
	if(ipop2.ne.0.and.irow.ne.nrow-1)vb2=dv_mod(ipop2)

!write(*,*)' ycur=',ycur,' ya2=',ya2,' yb2=',yb2
!write(*,*)' va2=',va2,' vb2=',vb2


	sa2b1=sqrt((za2-zb1)*(za2-zb1)+(x2-x1)*(x2-x1))
	sa1b2=sqrt((za1-zb2)*(za1-zb2)+(x2-x1)*(x2-x1))
	zt(1)=za1
	xt(1)=x1
	dvt(1)=va1
	zt(2)=zb2
	xt(2)=x2
	dvt(2)=vb2
	if(sa1b2.gt.sa2b1)then
		zt(1)=za2
		xt(1)=x2
		dvt(1)=va2
		zt(2)=zb1
		xt(2)=x1
		dvt(2)=vb1
	end if
	zper=zt(1)+(zt(1)-zt(2))/(xt(1)-xt(2))*(xcur-xt(1))
	if(zcur.ge.zper.and.sa1b2.le.sa2b1)then
		zt(3)=zb1
		xt(3)=x1
		dvt(3)=vb1
	else if(zcur.lt.zper.and.sa1b2.le.sa2b1)then
		zt(3)=za2
		xt(3)=x2
		dvt(3)=va2
	else if(zcur.ge.zper.and.sa1b2.gt.sa2b1)then
		zt(3)=zb2
		xt(3)=x2
		dvt(3)=vb2
	else if(zcur.lt.zper.and.sa1b2.gt.sa2b1)then
		zt(3)=za1
		xt(3)=x1
		dvt(3)=va1
	end if
	do i=1,3
		amatr(i,1)=xt(i)-xcur
		amatr(i,2)=zt(i)-zcur
		amatr(i,3)=1.
		bmatr(i)=dvt(i)
!write(*,*)(amatr(i,j),j=1,3),bmatr(i)
	end do
	call kram3(amatr,bmatr,sol)
	dvv(il)=sol(3)
!write(*,*)' dvv=',dvv(il)




end do

umn=wgt(1)+((wgt(2)-wgt(1))/(y2-y1))*(ycur-y1)
dv =dvv(1)+((dvv(2)-dvv(1))/(y2-y1))*(ycur-y1)

!write(*,*)' umn=',umn,' dv=',dv

return
end