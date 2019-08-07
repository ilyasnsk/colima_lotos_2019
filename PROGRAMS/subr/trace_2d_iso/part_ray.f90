subroutine part_ray(xzt,yzt,xst,yst,part1,part2)
real sray(1000)
common/ray/ nodes,xray(1000),yray(1000)
common/ray_part/ npart,xpart(1000),ypart(1000),spart(1000),kod_part(1000)

sss=0
sray(1)=0
do inode=1,nodes-1
	x1=xray(inode)
	x2=xray(inode+1)

	y1=yray(inode)
	y2=yray(inode+1)

	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))

	sss=sss+ds
	sray(inode+1)=sss
	!write(*,*)' t=',ttt,' s=',sss,' vvv=',vvv
end do

sA=part1*sss
sB=part2*sss
stotal=sss

n_A=0
dsmin=9999
do inode=1,nodes
	ds=abs(sray(inode)-sA)
	if(ds.gt.dsmin) cycle
	dsmin=ds
	n_A=inode
end do

n_B=0
dsmin=9999
do inode=1,nodes
	ds=abs(sray(inode)-sB)
	if(ds.gt.dsmin) cycle
	dsmin=ds
	n_B=inode
end do

npart=0
do inode=n_A,n_B

	x1=xray(inode)
	y1=yray(inode)
	s1=sray(inode)-sray(n_A)

	npart=npart+1
	xpart(npart)=x1
	ypart(npart)=y1
	spart(npart)=s1
	kod_part(npart)=inode
end do

!do i=1,npart
	!write(*,*)i,spart(i)/stotal,xpart(i),ypart(i) 
!end do

return
end