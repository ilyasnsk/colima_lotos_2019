subroutine remeshing(nrefl)
common/ray_param/ds_ini,ds_segm_min,bend_min0,bend_max0
common/ray/ nodes,xray(1000),yray(1000)
real xtmp(1000),ytmp(1000)

nmic0=200
dsmic0=ds_ini/nmic0

xtmp(1)=xray(1)
ytmp(1)=yray(1)

scur=0
ntmp=1
nrefl_new=0
do inode=1,nodes-1
	x1=xray(inode)
	x2=xray(inode+1)
	y1=yray(inode)
	y2=yray(inode+1)
	ds=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))
	nmic=int_best(ds/dsmic0)
	dsmic=ds/nmic
	dss=0
	do imic=1,nmic
		scur=scur+dsmic
		if(inode+1.eq.nrefl) then
			dss=dss+dsmic
			if(ds-dss.lt.ds_ini/2.) then
				ntmp=ntmp+1
				xtmp(ntmp)=x2
				ytmp(ntmp)=y2
				nrefl_new=ntmp
				yrefl=yyy_surf(x2)
				!write(*,*)' reflected:',ntmp,x2,y2,yrefl

				exit
			end if
		end if
		!write(*,*)imic,scur
		if(scur.ge.ds_ini) then
			ntmp=ntmp+1
			sss=dsmic*imic
			xtmp(ntmp)=x1+((x2-x1)/(ds))*sss
			ytmp(ntmp)=y1+((y2-y1)/(ds))*sss
			!write(*,*)' ntmp=',ntmp,' x=',xtmp(ntmp),' y=',ytmp(ntmp)
			scur=0
		end if
	end do

end do

ntmp=ntmp+1
xtmp(ntmp)=xray(nodes)
ytmp(ntmp)=yray(nodes)
!write(*,*)' ntmp=',ntmp

!open(23,file='ray_tmp.dat')
!do i=1,ntmp
	!write(23,*)xtmp(i),ytmp(i)
!end do
!close(23)


!write(*,*)' nodes=',nodes,' ntmp=',ntmp

nodes=ntmp
xray=xtmp
yray=ytmp
!write(*,*)' nodes=',nodes,' ntmp=',ntmp,' nrefl_new=',nrefl_new
nrefl=nrefl_new


return
end