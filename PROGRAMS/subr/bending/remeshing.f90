subroutine remeshing()

common/ray/ nodes,xray(1000),yray(1000),zray(1000)
real xtmp(1000),ytmp(1000),ztmp(1000)
common/ray_param/ds_ini


ntmp=nodes
xtmp=xray
ytmp=yray
ztmp=zray

nodes=0

do inode=1,ntmp-1
	x1=xtmp(inode)
	x2=xtmp(inode+1)

	y1=ytmp(inode)
	y2=ytmp(inode+1)
	
	z1=ztmp(inode) 
	z2=ztmp(inode+1) 

	dz1=ztmp(inode) 
	dz2=ztmp(inode+1) 

	nodes=nodes+1
	xray(nodes)=x1
	yray(nodes)=y1
	zray(nodes)=z1

	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
	!write(*,*)' normal: ',ds,xray(nodes),yray(nodes),zray(nodes)

	if(ds.gt.ds_ini*1.5) then
		nmid = ds/ds_ini + 1
		dsmid = ds/nmid
		do imid=1,nmid-1

			xmid=x1+((x2-x1)/ds)*(dsmid*imid)
			ymid=y1+((y2-y1)/ds)*(dsmid*imid)
			zmid=z1+((z2-z1)/ds)*(dsmid*imid)

			nodes=nodes+1
			xray(nodes)=xmid
			yray(nodes)=ymid
			zray(nodes)=zmid
		!write(*,*)'  mid: ',dsmid,xray(nodes),yray(nodes),zray(nodes)
		end do
	end if
end do

nodes=nodes+1
xray(nodes)=xtmp(ntmp)
yray(nodes)=ytmp(ntmp)
zray(nodes)=ztmp(ntmp)


return
end