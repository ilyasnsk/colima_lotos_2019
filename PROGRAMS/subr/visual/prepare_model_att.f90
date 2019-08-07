subroutine prepare_model_att(ar,md,igr,ips)

character*1 gr,ps
character*8 ar,md,line

common/pi/pi,per
common/center/fi0,tet0
common/grid/nornt,ornt(4),sinal,cosal,&
		nlev,ylev(300),ntop(300),nobr,&
		xtop(20000,300), ztop(20000,300),n_pop(20000,300),&
		att_mod(100000)

write(gr,'(i1)')igr
write(ps,'(i1)')ips


!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=573)line
	if(line.eq.'ORIENTAT') goto 574
end do
573 continue
write(*,*)' cannot find ORIENTATIONS in MAJOR_PARAM.DAT!!!'
pause
574 read(1,*)nornt
read(1,*)(ornt(i),i=1,nornt)
close(1)
!******************************************************************

!write(*,*)' fi0=',fi0,' tet0=',tet0
orient=ornt(igr)
sinal=sin(orient*per)
cosal=cos(orient*per)
!write(*,*)' orient=',orient


open(1,file='../../../DATA/'//ar//'/'//md//'/data/lev_att'//gr//ps//'.dat')
i=0
722 i=i+1
	read(1,*,end=721)n,ylev(i)
	!write(*,*)i,n,ylev(i)
	goto 722
721 nlev=i-1
close(1)



ntop=0
open(1,file='../../../DATA/'//ar//'/'//md//'/data/gr_att'//gr//ps//'.dat')
do n=1,nlev
	read(1,*,end=556)ntop(n)
	!write(*,*)' y=',ylev(n),' ntop=',ntop(n)
	do i=1,ntop(n)
		read(1,*)xtop(i,n),ztop(i,n),n_pop(i,n)
	end do
end do
556		close(1)

open(1,file='../../../DATA/'//ar//'/'//md//'/data/obr_att'//gr//ps//'.dat')
read(1,*) nobr
if(nobr.gt.100000) then
	write(*,*)' nobr=',nobr
	write(*,*)' One should increase size of dv_mod'
	pause
end if
close(1)
!write(*,*)' nobr=',nobr

open(1,file='../../../DATA/'//ar//'/'//md//'/data/att_result_'//gr//ps//'.dat')
do i=1,nobr
	read(1,*)  att_mod(i)
end do
close(1)

return
end

