subroutine prep_board_dv(ar,md)
character*8 ar,md
common/brd_dv/xbr1(2),xbr2(2),dxbr1(2),dxbr2(2),&
			  ybr1(2),ybr2(2),dybr1(2),dybr2(2),&
			  zbr1(2),zbr2(2),dzbr1(2),dzbr2(2),amp(4)

open(1,file='../../../DATA/'//ar//'/'//md//'/anomaly.dat')
read(1,*)				
read(1,*)				
read(1,*)amp(1),amp(3)				
read(1,*)xbr1(1),xbr2(1),dxbr1(1),dxbr2(1)				
read(1,*)ybr1(1),ybr2(1),dybr1(1),dybr2(1)				
read(1,*)zbr1(1),zbr2(1),dzbr1(1),dzbr2(1)
read(1,*)amp(2),amp(4)				
read(1,*)xbr1(2),xbr2(2),dxbr1(2),dxbr2(2)				
read(1,*)ybr1(2),ybr2(2),dybr1(2),dybr2(2)				
read(1,*)zbr1(2),zbr2(2),dzbr1(2),dzbr2(2)				
close(1)

write(*,*)' amp1=',amp(1),' amp2=',amp(2)


return
end