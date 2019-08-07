character*8 ar,md,line
character*1 it
integer k_bad(2,50000)


open(1,file='../../../model.dat')
read(1,'(a8)')ar		! code of the area
read(1,'(a8)')md		! code of the area
read(1,*)niter		! code of the area
close(1)

!******************************************************************
key_ft1_xy2=1
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
    read(1,'(a8)',end=513)line
    if(line.eq.'GENERAL ') goto 514
end do
513 continue
write(*,*)' cannot find GENERAL INFORMATION in MAJOR_PARAM.DAT!!!'
pause
514 continue
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*,end=441,err=441)key_ft1_xy2
read(1,*,end=441,err=441)key_true1
441 close(1)
!******************************************************************


open(1,file='../../../PROGRAMS/3_VISUAL/variance_reduction/SET.dat')
read(1,*)norm		! norm for dispersion
read(1,*)alim_p
read(1,*)alim_s
close(1)
write(*,*)' ar=',ar,' md=',md,' niter=',niter

k_bad=0
nbad=0

open(21,file='../../../PICS/'//ar//'/'//md//'/info_resid.dat')

do iter=1,niter
    write(it,'(i1)')iter
    nzt=0
    nray=0

    dtot_p=0
    dtot_s=0

    ntot_p=0
    ntot_s=0

    err_loc=0


	open(1,file='../../../DATA/'//ar//'/'//md//'/data/rays'//it//'.dat')
	open(11,file='../../../PICS/'//ar//'/'//md//'/final_resid.dat')

	! Read the sources:
	872	continue
	        if(key_true1.eq.1) read(1,*,end=871)xtrue,ytrue,ztrue
               read(1,*,end=871)xzt,yzt,zzt,nkrat

                if(key_true1.eq.1) then
                    err = sqrt ( (xtrue-xzt)**2 + (ytrue-yzt)**2 + (ztrue-zzt)**2 )
                    err_loc = err_loc + err
                end if


		nzt=nzt+1
		write(11,*)nzt,xzt,yzt,zzt,nkrat
		nray=nray+nkrat
		ddd_p=0
		nnn_p=0
		ddd_s=0
		nnn_s=0
		do ikr=1,nkrat
			read(1,*)ips,ist,tobs,tmod
			dt=tobs-tmod
			write(11,'(2i4,f8.3)')ips,ist,dt
			if(iter.eq.1) then
				alim=alim_p
				if(ips.eq.2) alim=alim_s
				if(abs(dt).gt.alim) then
					nbad=nbad+1
					k_bad(1,nbad)=nzt
					k_bad(2,nbad)=ikr
					cycle
				end if
			else
				do i=1,nbad
					if(k_bad(1,i).eq.nzt.and.k_bad(2,i).eq.ikr) goto 874
				end do
			end if
			!write(*,*)' ips=',ips,' dt=',dt

			if(ips.eq.1) then
				if(norm.eq.2) then
					ddd_p = ddd_p + dt*dt
				else if(norm.eq.1) then
					ddd_p = ddd_p + abs(dt)
				end if
				nnn_p = nnn_p + 1
			else
				if(norm.eq.2) then
					ddd_s = ddd_s + dt*dt
				else if(norm.eq.1) then
					ddd_s = ddd_s + abs(dt)
				end if
				nnn_s = nnn_s + 1
			end if
874			continue
		end do
		dcur_p=0
		dcur_s=0
		if(norm.eq.2) then
			dcur_p=sqrt(ddd_p/nnn_p)
			if(nnn_s.ne.0) dcur_s=sqrt(ddd_s/nnn_s)
		else if(norm.eq.1) then
			dcur_p=ddd_p/nnn_p
			if(nnn_s.ne.0) dcur_s=ddd_s/nnn_s
		end if

		write(11,*)' dcur_p=',dcur_p,' dcur_s=',dcur_s

		!write(*,*)' dcur_p=',dcur_p,' nnn_p=',nnn_p
		!write(*,*)' dcur_s=',dcur_s,' nnn_s=',nnn_s

		
		dtot_p=dtot_p+ddd_p
		ntot_p=ntot_p+nnn_p
		dtot_s=dtot_s+ddd_s
		ntot_s=ntot_s+nnn_s
		goto 872
871	close(1)
	close(11)
	if(norm.eq.2) then
		dtot_p=sqrt(dtot_p/ntot_p)
		dtot_s=sqrt(dtot_s/ntot_s)
	else if(norm.eq.1) then
		dtot_p=(dtot_p/ntot_p)
		dtot_s=(dtot_s/ntot_s)
	end if
	if(iter.eq.1) then
		write(*,*)' nbad=',nbad
		dtot1_p=dtot_p
		dtot1_s=dtot_s
	end if
	red_p=100*(dtot1_p-dtot_p)/dtot1_p
	red_s=100*(dtot1_s-dtot_s)/dtot1_s
	write(*,*)' iter=',iter,' dtot_p=',dtot_p,' red=',red_p
	write(*,*)' iter=',iter,' dtot_s=',dtot_s,' red=',red_s
        if(key_true1.eq.1)write(*,*)' Source mislocation:',err_loc/nzt
	write(*,*)'___________________________________________________'
	write(21,*)' iter=',iter,' dtot_p=',dtot_p,' red=',red_p
	write(21,*)' iter=',iter,' dtot_s=',dtot_s,' red=',red_s
        if(key_true1.eq.1)write(21,*)' Source mislocation:',err_loc/nzt
	write(21,*)'___________________________________________________'


end do
write(*,*)' nsrces=',nzt,' nray_p=',ntot_p,' nray_s=',ntot_s
write(21,*)' nsrces=',nzt,' nray_p=',ntot_p,' nray_s=',ntot_s

close(21)
close(11)
stop
end