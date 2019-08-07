function anomaly(xxx,yyy,zzz,ips)
common/kod_anom/n_anomaly

dv=0

!write(*,*)' n_anomaly=',n_anomaly

if(n_anomaly.eq.1)then
	dv = dv_board(xxx,yyy,zzz,ips)
else if(n_anomaly.eq.2)then
	dv = hor_anom(xxx,yyy,zzz,ips)
else if(n_anomaly.eq.3)then
	dv=vert_anom(xxx,yyy,zzz,ips)
else if(n_anomaly.eq.4)then
	dv=vert_brd(xxx,yyy,zzz,ips)
!else if(n_anomaly.eq.5)then
!	dv=vert_PS_separ(xxx,yyyy,zzz,ips)
end if

anomaly=dv

return
end