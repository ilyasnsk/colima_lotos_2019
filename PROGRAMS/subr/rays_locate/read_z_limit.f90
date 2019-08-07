subroutine read_z_lim(ar,md)
character*8 md,ar,line
common/zlimit/zbase,nar,npar(5),zmar(5),flim(5,100),tlim(5,100)

nar=0
zbase=100000

! Maximal depth is defined as deepest level in the reference table:
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,*,end=553)line
	if(line.eq.'REF_PARA') goto 554
end do
553 continue
write(*,*)' cannot find REF_PARAM in MAJOR_PARAM.DAT!!!'
goto 555

554 read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)nstep		
do i=1,nstep
	read(1,*)
end do
read(1,*)zztmax
		

zbase=zztmax
555 continue
close(1)

write(*,*)' ar=',ar,' md=',md

! Maximal depth can be also defined in polylines
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
    read(1,*,end=563)line
    !write(*,*)line
    if(line.eq.'SRCE_LIM') goto 564
end do
563 continue
!write(*,*)' cannot find SOURCE_LIMITS in MAJOR_PARAM.DAT!!!'
goto 565

564 continue
zbase2=100000
nar=0
read(1,*,err=1)zbase2
!write(*,*)' zbase2=',zbase2
2 read(1,*,err=1)np,zmax
!write(*,*)' np=',np,' zmax=',zmax
nar=nar+1
npar(nar)=np
zmar(nar)=zmax
do i=1,np
	read(1,*) flim(nar,i), tlim(nar,i)
end do
goto 2
1 continue
!write(*,*)' nar=',nar

if(zbase2.lt.zbase) zbase=zbase2

565 continue
close(1)

write(*,*)' zbase=',zbase
!pause
return
end