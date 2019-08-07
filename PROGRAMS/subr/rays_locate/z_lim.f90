function z_lim(fi,tet)
common/zlimit/zbase,nar,npar(5),zmar(5),flim(5,100),tlim(5,100)

z_lim=zbase

!write(*,*)' zbase=',zbase,' nar=',!
if(nar.eq.0) return

do ian=1,nar
	zmax=zmar(ian)
	!write(*,*)'zmar=',zmar(ian),'npar=',npar(ian)

	icrup=0
	do inod=2,npar(ian)
		!write(*,*)' fi=',anom(ian,inod,1),' tet=',anom(ian,inod,2)
		f1=flim(ian,inod-1)
		f2=flim(ian,inod)
		t1=tlim(ian,inod-1)
		t2=tlim(ian,inod)
		if(f1.eq.f2) cycle
		if ((fi-f1)*(fi-f2).gt.0.) cycle
		!write(*,*)' f1=',f1,' f2=',f2
		!write(*,*)' t1=',t1,' t2=',t2
		if ((fi-f1).eq.0.) then
			in=inod
667			continue
			if(in.gt.2) then
				f0=flim(ian,in-2)
				if(f0.eq.f1) then
					in=in-1
					goto 667
				end if
			else
				f0=flim(ian,npar(ian)-1)
			end if
			!write(*,*)' i=',inod,' f0=',f0,' f1=',f1,' f2=',f2
			!pause
			if ((f0-f1)*(f2-f1).gt.0.) cycle
		end if
		if ((fi-f2).eq.0.) cycle
		tt=t1+((t2-t1)/(f2-f1))*(fi-f1)
		if (tt.lt.tet) cycle
		!write(*,*)' tt=',tt
		icrup=icrup+1
	end do
	icrup2=int(icrup/2)*2
	!write(*,*)' ian=',ian,' icrup=',icrup,' icrup2=',icrup2
	if(icrup.eq.icrup2) cycle
	z_lim=zmar(ian)
	exit
end do

return
end