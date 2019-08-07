subroutine refmod_all(distkm,depzt,ips, nall,tall,hall,aall)
real time2(30),alfa2(30),hmax2(30)
real time3(2,30),alfa3(2,30),hmax3(2,30)
integer ikr3(2),k_test(2,30)
real tall(20),hall(20),aall(20)

common/reftable/izttab,ntab(2,200),hzttab(2,200),ttab(2,200,5000),rtab(2,200,5000),etab(2,200,5000),atab(2,200,5000)
common/deriv/dtdd,dtdz
common/crust_par/avmoho,vcr1_p,vcr2_p,vcr1_s,vcr2_s,vsrf_p,vsrf_s,acrs,acrp


REAL PI/3.1415926/
PER=PI/180.

nall=0
time=-900
depth=-900
alfa=-900

!write(*,*)' ips=',ips,' izttab=',izttab
do iz=1,izttab
	
	z1=hzttab(ips,iz)
	z2=hzttab(ips,iz+1)
	if((depzt-z1)*(depzt-z2).le.0.) goto 1
end do
write(*,*)' cannot find for the following z:'
write(*,*)' zzt=',depzt
pause
return
1 continue
!write(*,*)' depzt=',depzt,' z1=',z1,' z2=',z2

alfa2=-999
time2=-999
hmax2=-999

istep=0
3 istep=istep+1
	izz=iz
	if(istep.eq.2)izz=iz+1
	nn1=1
	nn2=ntab(ips,izz)
		!write(*,*)' nn1=',nn1,' nn2=',nn2

	ikrat=0
	do iray=nn1,nn2-1
		d1=etab(ips,izz,iray)
		d2=etab(ips,izz,iray+1)
		!write(*,*)' d1=',d1,' d2=',d2,' distkm=',distkm
		if((distkm-d1)*(distkm-d2).gt.0.) cycle
		a1=atab(ips,izz,iray)
		a2=atab(ips,izz,iray+1)
		t1=ttab(ips,izz,iray)
		t2=ttab(ips,izz,iray+1)
        !write(*,*)' iray=',iray
		!write(*,*)' FOUND! d1=',d1,' d2=',d2,' distkm=',distkm
		!write(*,*)' t1=',t1,' t2=',t2,ikrat
		h1=rtab(ips,izz,iray)
		h2=rtab(ips,izz,iray+1)
        if(ikrat+1.gt.20) exit
		ikrat=ikrat+1
		alfa3(istep,ikrat)=a1+((a2-a1)/(d2-d1))*(distkm-d1)
		time3(istep,ikrat)=t1+((t2-t1)/(d2-d1))*(distkm-d1)
		hmax3(istep,ikrat)=h1+((h2-h1)/(d2-d1))*(distkm-d1)
		!write(*,*)' distkm=',distkm,' d1=',d1,' d2=',d2
!write(*,*)time3(istep,ikrat),alfa3(istep,ikrat),hmax3(istep,ikrat)
		!pause
	end do
	ikr3(istep)=ikrat
21	continue

if(istep.eq.1) goto 3

if(ikr3(1).eq.0.or.ikr3(2).eq.0) return

ikmin=ikr3(1)
if(ikr3(2).lt.ikr3(1)) ikmin=ikr3(2)

k_test=0

ik3=0

11 continue
hmin=99999
do ik1=1,ikr3(1)
	do ik2=1,ikr3(2)
		if(k_test(1,ik1).eq.1.and.k_test(1,ik1).eq.1) cycle
		dh=abs(hmax3(1,ik1)-hmax3(2,ik2))
		if(dh.gt.hmin) cycle
		i1=ik1
		i2=ik2
		hmin=dh
	end do
end do

ik3=ik3+1

a1=alfa3(1,i1)
a2=alfa3(2,i2)
t1=time3(1,i1)
t2=time3(2,i2)
h1=hmax3(1,i1)
h2=hmax3(2,i2)
!write(*,*)' i1=',i1,' i2=',i2
!write(*,*)' h1=',h1,' h2=',h2
!write(*,*)' t1=',t1,' t2=',t2
!write(*,*)' a1=',a1,' a2=',a2

k_test(1,i1)=1
k_test(2,i2)=1

aall(ik3)=a1+((a2-a1)/(z2-z1))*(depzt-z1)
tall(ik3)=t1+((t2-t1)/(z2-z1))*(depzt-z1)
hall(ik3)=h1+((h2-h1)/(z2-z1))*(depzt-z1)

!write(*,*)' ik3=',ik3,' a=',aall(ik3),' t=',tall(ik3),' h=',hall(ik3)
!pause

if(ik3.lt.ikmin) goto 11

nall=ik3

return
end