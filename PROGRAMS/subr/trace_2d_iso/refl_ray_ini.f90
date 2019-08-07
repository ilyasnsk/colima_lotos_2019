subroutine refl_ray_ini(xzt,yzt,xst,yst,ips, nrefl,tout)

real xtmp(1000),ytmp(1000)
common/ray/ nodes,xray(1000),yray(1000)
common/ray_param/ds_ini,ds_part_min,bend_min0,bend_max0

yzt2 = yyy_surf(xzt)*2 - yzt
yst2 = yyy_surf(xst)*2 - yst

write(*,*)' yzt2=',yzt2,' yst2=',yst2

aaa1=(yzt-yst2)/(xzt-xst)
bbb1=(xzt*yst2-xst*yzt)/(xzt-xst)

aaa2=(yzt2-yst)/(xzt-xst)
bbb2=(xzt*yst-xst*yzt2)/(xzt-xst)

xxx = (bbb2-bbb1) / (aaa1-aaa2)
yyy1 = aaa1 * xxx + bbb1
yyy2 = aaa2 * xxx + bbb2

write(*,*)' xxx=',xxx,' yyy=',yyy1,yyy2

xrefl=xxx
yrefl=yyy_surf(xxx)

write(*,*)' xrefl=',xrefl,' yrefl=',yrefl


call streight_line(xzt,yzt,xrefl,yrefl,ips, tout1)

nrefl=nodes

ntmp=nodes
xtmp=xray
ytmp=yray

call streight_line(xrefl,yrefl,xst,yst,ips, tout2)

do i=2,nodes
	ntmp=ntmp+1
	xtmp(ntmp)=xray(i)
	ytmp(ntmp)=yray(i)
end do

nodes=ntmp
xray=xtmp
yray=ytmp

tout=tout1+tout2

!do i=1,nodes
	!write(*,*)i,xray(i),yray(i)
!end do

return
end