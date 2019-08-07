subroutine prepare_noise()

USE DFPORT

real hist(2,100),ahist(1000)
character*1 it_ini
integer*2 k12,k22,k32,k42

common/histo/istep,xconv(2000),xdir1,xdir2,dxdir
common/noise_1/a_rand(2),xmid(2)


!******************************************************************
kod_ini1_mod2=1
nhistog=60
nsample=10000


call gettim(k12,k22,k32,k42)
k4=k42
CALL SRAND(k4)


open(1,file='SET.DAT')
read(1,*)a_rand_p, a_rand_s
read(1,*)xmid(1), xmid(2)
close(1)


open(1,file='../../../COMMON/noise_histo.dat',status='old',err=221)
! Read histogram of noise distribution:
read(1,*)nconvert
read(1,*)ds
do i=1,100
	read(1,*,end=1)hist(1,i),hist(2,i)
end do
pause ' Warning: nhist > 100'
1 nhist=i-1
close(1)
goto 222
221 write(*,*)' file COMMON/noise_histo.dat does not exist'
write(*,*)' create this file or use the option "2" for the noise type'
stop
222 continue

write(*,*)' nconvert=',nconvert,' nhist=',nhist
write(*,*)' Noise: amp=',a_rand(1),a_rand(2)

amax=0
do i=1,nhist
	if(hist(2,i).gt.amax) amax=hist(2,i)
end do


xdir1=0
xdir2=1
dxdir0=(xdir2-xdir1)/nconvert

! Integration

sum=0
do i=1,nhist-1
	x1=hist(1,i)
	y1=hist(2,i)
	x2=hist(1,i+1)
	y2=hist(2,i+1)
	sum=sum+(x2-x1)*(y1+y2)/2
end do

!write(*,*)' sum=',sum

porog=sum/nconvert

istep=1


xconv(istep)=hist(1,1)

dsum=0
do i=1,nhist-1
	x1=hist(1,i)
	x2=hist(1,i+1)
	y1=hist(2,i)
	y2=hist(2,i+1)
	nstep=(x2-x1)/ds
	dss=(x2-x1)/nstep
	dydx=(y2-y1)/(x2-x1)
	x3=x1
	y3=y1
	do ii=1,nstep
		x4=x3+dss
		y4=y1+dydx*(x4-x1)
		dsum=dsum+(x4-x3)*(y4+y3)/2
		x3=x4
		y3=y4

		if(dsum.lt.porog) cycle

		!write(*,*)' x4=',x4,' dsum=',dsum,' porog=',porog
		istep=istep+1
		xconv(istep)=x4
		dsum=0
	end do
end do
istep=istep+1
xconv(istep)=hist(1,nhist)

dxdir=(xdir2-xdir1)/(istep-1)


a_rand(1)=1

disp=0
do i=1,nsample
	a=our_noise(1)
	disp=disp+abs(a)
end do
disp=disp/nsample
!write(*,*)' disp P=',disp

a_rand(1)=a_rand_p/disp

amin=-a_rand_p*4
amax=a_rand_p*4
dhistog=(amax-amin)/nhistog

disp=0
ahist=0
do i=1,nsample
    a=our_noise(1)
    disp=disp+abs(a)
    if((a-amin)*(a-amax).gt.0) cycle
    do ih=1,nhistog
        a1=amin+(ih-1)*dhistog
        a2=amin+(ih)*dhistog
        if((a-a1)*(a-a2).le.0.)exit
    end do
    ahist(ih)=ahist(ih)+1
end do
disp=disp/nsample
write(*,*)' disp P=',disp

open(11,file='data_out/hist_p.bln')
write(11,*)nhistog
do ih=1,nhistog
    ah=amin+(ih-0.5)*dhistog
    values=ahist(ih)*100/nsample
    !write(*,*)ih,ah,ahist(ih)
    write(11,*)ah,values
enddo
close(11)


!*****************************************************

a_rand(2)=1
disp=0
do i=1,nsample
    a=our_noise(2)
    disp=disp+abs(a)
end do
disp=disp/nsample
!write(*,*)' disp S=',disp

a_rand(2)=a_rand_s/disp

amin=-a_rand_s*4
amax=a_rand_s*4
dhistog=(amax-amin)/nhistog
!write(*,*)' amin=',amin,' amax=',amax,' nhistog=',nhistog,' dhistog=',dhistog

disp=0
ahist=0
do i=1,nsample
    a=our_noise(2)
    disp=disp+abs(a)
    if((a-amin)*(a-amax).gt.0) cycle
    do ih=1,nhistog
        a1=amin+(ih-1)*dhistog
        a2=amin+(ih)*dhistog
        if((a-a1)*(a-a2).le.0.)exit
    end do
    ahist(ih)=ahist(ih)+1
    !if(ih.eq.41) write(*,*)a
end do


open(11,file='data_out/hist_s.bln')
write(11,*)nhistog
do ih=1,nhistog
    ah=amin+(ih-0.5)*dhistog
    values=ahist(ih)*100/nsample
    !write(*,*)ih,ah,ahist(ih)
    write(11,*)ah,values
enddo
close(11)
disp=disp/nsample
write(*,*)' disp S=',disp



return
end