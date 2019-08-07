subroutine prepare_ref(ar,md)
character*8 ar,md

common/reftable/izttab,ntab(2,200),hzttab(2,200),ttab(2,200,5000),rtab(2,200,5000),etab(2,200,5000),atab(2,200,5000)
common/zlimit/zbase,nar,npar(5),zmar(5),flim(5,100),tlim(5,100)


open(1,file='../../../DATA/'//ar//'/'//md//'/data/table.dat',form='binary')
izt=0
nrmax=0
34	continue
	izt=izt+1
	do ips=1,2
		read(1,end=35)hzttab(ips,izt),ntab(ips,izt)
		if(ntab(ips,izt).gt.nrmax) nrmax=ntab(ips,izt)
		!write(*,*)' i=',izt,' z=',hzttab(ips,izt),' ips=',ips,' nref=',ntab(ips,izt)
		do i=1,ntab(ips,izt)
			read(1)etab(ips,izt,i),ttab(ips,izt,i),atab(ips,izt,i),rtab(ips,izt,i)
			!write(*,*)etab(ips,izt,i),atab(ips,izt,i)
		end do
		!write(*,*)izt,' z=',hzttab(ips,izt),' n=',ntab(ips,izt),' d=',etab(ips,izt,ntab(ips,izt))
	end do
goto 34
35 close(11)
izttab=izt-1
if(zbase.gt.hzttab(1,izttab))zbase=hzttab(1,izttab)
!write(*,*)' izttab=',izttab,' zbase=',zbase
return 
end