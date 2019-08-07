subroutine read_anom(ar,md)
character*8 ar,md
common/kod_anom/n_anomaly

open(1,file='../../../DATA/'//ar//'/'//md//'/anomaly.dat')
read(1,*) n_anomaly
close(1)
write(*,*)' ar=',ar,'   kod of anom.=',n_anomaly


if(n_anomaly.eq.1)then
    write(*,*)' n_anomaly=1, : CHECKERBOARD'
	call prep_board_dv(ar,md)
else if(n_anomaly.eq.2)then
    write(*,*)' n_anomaly=2, : FREE HORIZONTAL ANOMALIES'
	call read_hor_an (ar,md)
else if(n_anomaly.eq.3)then
    write(*,*)' n_anomaly=3, : FREE VERTICAL ANOMALIES'
	call read_vert_an(ar,md)
else if(n_anomaly.eq.4)then
    write(*,*)' n_anomaly=4, : VERTICAL CHECKERBOARD'
	call read_vert_brd(ar,md)
!else if(n_anomaly.eq.5)then
!	call read_PS_separ(ar,md)
end if


return
end