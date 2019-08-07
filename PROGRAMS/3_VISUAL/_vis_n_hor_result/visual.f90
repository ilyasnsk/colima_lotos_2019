USE DFPORT

! HORIZONTAL !!!
! NODES !!!!!!!

character*4 dsaa/'DSAA'/
character*8 ar,md,line
character*2 lv
character*1 ps ,rm,it
character*20 scale_line, scale_line2
allocatable dvan(:,:),vvv(:,:),v1tmp(:,:),v2tmp(:,:),vabs(:,:,:)
real hlev(20), vaver(2,20)
integer nrps(2)
real fzzt(10000,100),tzzt(10000,100),zzzt(10000,100)
integer nzzt(100),line1_rgb(3),line2_rgb(3),kdot_rgb(3)
integer*2 izzz

common/pi/pi,per
common/center/fi0,tet0
common/keys/key_ft1_xy2

one=1.e0
pi=asin(one)*2.e0
per=pi/180.e0
rz=6371.

w_limit=0.2

igr=1

open(1,file='../../../model.dat')
read(1,'(a8)')ar
read(1,'(a8)')md
read(1,*)iter		
close(1)
write(it,'(i1)')iter

write(*,*)
write(*,*)' ***********************************************'
write(*,*)' VISUALISATION in horizontal sections: '
write(*,*)' ar=',ar,' md=',md,' iter=',iter


i=system('mkdir ..\..\..\TMP_files\hor')


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
441 close(1)
!******************************************************************

key_preview=0
open(1,file='../../../preview_key.txt')
read(1,*,end=771)key_preview
771 close(1)

if(key_preview.ne.0) then
	i=system('mkdir ..\..\..\PICS\'//ar//'\'//md//'\IT'//it)
	open(1,file='../../../DATA/'//ar//'/config.txt')
	read(1,*) 
	read(1,*) npix_x0,npix_y0
	read(1,*)tick_x,tick_y
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*) 
	read(1,*)scale_line
	read(1,*)amp_min,amp_max
	read(1,*)scale_line2
	read(1,*)vpvs_min,vpvs_max
	close(1)

	i=system('copy ..\..\..\COMMON\visual_exe\visual.exe layers.exe')
	i=system('copy ..\..\..\COMMON\scales_scl\'//scale_line//' '//scale_line)
	i=system('copy ..\..\..\COMMON\scales_scl\'//scale_line2//' '//scale_line2)

end if

!******************************************************************
open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
do i=1,10000
	read(1,'(a8)',end=573)line
	if(line.eq.'ORIENTAT') goto 574
end do
573 continue
write(*,*)' cannot find ORIENTATIONS in MAJOR_PARAM.DAT!!!'
pause
574 read(1,*)nornt
close(1)
!******************************************************************

ngr1=1
ngr2=nornt


kod_av_bias=0
kod_apriori=0
ind_srce=1


write(*,*)' ar=',ar,' md=',md
write(it,'(i1)')iter

!******************************************************************
if(key_ft1_xy2.eq.1) then
    open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
    do i=1,10000
        read(1,'(a8)',end=553)line
        if(line.eq.'AREA_CEN') goto 554
    end do
    553 continue
    write(*,*)' cannot find AREA CENTER in MAJOR_PARAM.DAT!!!'
    pause
    554 read(1,*)fi0,tet0
    close(1)
else
    fi0=0
    tet0=0
end if
!******************************************************************

open(2,file='../../../DATA/'//ar//'/sethor.dat')
read(2,*) nlev  
read(2,*) (hlev(i),i=1,nlev)  
read(2,*) fmap1,fmap2,dfmap,tmap1,tmap2,dtmap  
read(2,*) smaxx
read(2,*) ismth
close(2)

cost=cos(tet0*per)
if(npix_x0.ne.0) then
    npix_x=npix_x0
else
    npix_x=int(npix_y0*((fmap2-fmap1)/(tmap2-tmap1))*cost)
    npix_y=npix_y0
end if

if(npix_y0.ne.0) then
    npix_y=npix_y0
else
    npix_y=int(npix_x0*((tmap2-tmap1)/(fmap2-fmap1))/cost)
    npix_x=npix_x0
end if

write(*,*)' npix_x=',npix_x,' npix_y=',npix_y


if(kod_apriori.eq.1) then
	call read_ini_model(ar,md)
end if

call read_topo(ar)

rsmth=ismth+0.5
nfmap=int_best((fmap2-fmap1)/dfmap+1.)
ntmap=int_best((tmap2-tmap1)/dtmap+1.)
write(*,*)' nfmap=',nfmap,' ntmap=',ntmap
allocate(dvan(nfmap,ntmap),vvv(nfmap,ntmap),v1tmp(nfmap,ntmap),v2tmp(nfmap,ntmap))
allocate(vabs(2,nfmap,ntmap))


!call read_3D_mod_v(ar,iter-1)

call read_vref(ar,md)


open(1,file='../../../DATA/'//ar//'/'//md//'/data/numray1.dat')
read(1,*) nrps(1),nrps(2)
close(1)


if(ind_srce.ne.0) then
	nzzt=0
	nzt=0
       open(1,file='../../../DATA/'//ar//'/'//md//'/data/srces'//it//'.dat')
	872	read(1,*,end=871)fzt,tzt,zzt
	!call decsf(xzt,yzt,0.,fi0,tet0,fzt,tzt,h)
	nzt=nzt+1
	do ilev=1,nlev
		if(ilev.eq.1) then
			z1=-10
			z2=(hlev(1)+hlev(2))/2
		else if(ilev.eq.nlev) then
			z1=(hlev(nlev-1)+hlev(nlev))/2
			z2=hlev(nlev) + (hlev(nlev)-hlev(nlev-1))/2
		else 
			z1=(hlev(ilev-1)+hlev(ilev))/2
			z2=(hlev(ilev+1)+hlev(ilev))/2
		end if
		if((zzt-z1)*(zzt-z2).le.0) goto 995
	end do
	goto 872

995 continue
	!write(*,*)' zzt=',zzt,' ilev=',ilev,' z1=',z1,' z2=',z2
	nzzt(ilev)=nzzt(ilev)+1
	fzzt(nzzt(ilev),ilev)=fzt
	tzzt(nzzt(ilev),ilev)=tzt
	zzzt(nzzt(ilev),ilev)=zzt

	goto 872
871 close(1)
	DO ilev=1,nlev
		write(lv,'(i2)')ilev
		open(11,file='../../../TMP_files/hor/ztr'//lv//'.dat')
		write(*,*)' ilev=',ilev,' nzzt=',nzzt(ilev)
		do izzt=1,nzzt(ilev)
			write(11,*)fzzt(izzt,ilev),tzzt(izzt,ilev),zzzt(izzt,ilev)
		end do
		close(11)
	end do
end if


vaver=0
DO ilev=1,nlev
	zzz=hlev(ilev)
	write(lv,'(i2)')ilev
	write(*,*)' ilev=',ilev,' zzz=',zzz
	do ips=1,2
		v0=vrefmod(zzz,ips)
		if(nrps(ips).eq.0) cycle
		write(ps,'(i1)')ips
		dvan=0
		vvv=0
		do igr=ngr1,ngr2
			call prepare_model_v(ar,md,ips,iter,igr)

			do itet=1,ntmap
				ttt=(itet-1)*dtmap+tmap1+tet0
				!write(*,*)' itet=',itet,' ttt=',ttt
				!ttt=-7.45
				do ifi=1,nfmap
					fff=(ifi-1)*dfmap+fmap1+fi0
					!fff=fi0+1
					!fff=110.87
					dv=0
					www=0
					call dv_1_grid_v(fff,ttt,zzz,smaxx, dv,www)
					zlim_up=h_lim(fff,ttt)
					if (zzz.lt.zlim_up) www=0
					!write(*,*)' dv=',dv,' www=',www
					dvproc=100*dv/v0
					dvan(ifi,itet)=dvan(ifi,itet)+dvproc*www
					vvv(ifi,itet)=vvv(ifi,itet)+www
					!if(itet.eq.101) write(*,*)dv,www
				end do
			end do
		end do

        nonzer=0
        do ifi=1,nfmap
            do itet=1,ntmap
                vanm=-999.
                vabs(ips,ifi,itet)=-999
                if (vvv(ifi,itet).gt.w_limit) then
                    vanm=dvan(ifi,itet)/vvv(ifi,itet)
                    vabs(ips,ifi,itet)=v0*(1+0.01*vanm)
                    nonzer=nonzer+1
                    vaver(ips,ilev)=vaver(ips,ilev)+vabs(ips,ifi,itet)
                end if
                v1tmp(ifi,itet)=vanm
            end do
        end do
        vaver(ips,ilev)=vaver(ips,ilev)/nonzer
        write(*,*)' ips=',ips,' vaver=',vaver(ips,ilev),' nonzer=',nonzer


		do ifi=1,nfmap
			do itet=1,ntmap
				if(vvv(ifi,itet).lt.w_limit) cycle
				vanm=0.
				iv=0
				do iff=-ismth,ismth
					if (ifi+iff.lt.1) cycle
					if (ifi+iff.gt.nfmap) cycle
					do itt=-ismth,ismth
						if (itet+itt.lt.1) cycle
						if (itet+itt.gt.ntmap) cycle
						if(vvv(ifi+iff,itet+itt).lt.w_limit) cycle
						rr=iff*iff+itt*itt
						r=sqrt(rr)
						if(r.gt.rsmth) cycle
						iv=iv+1
						vanm=vanm+v1tmp(ifi+iff,itet+itt)
					end do
				end do
				v2tmp(ifi,itet)=vanm/iv
			end do
		end do

		aver=0
		naver=0
		do ifi=1,nfmap
			do itet=1,ntmap
				vanom=-999
				if(vvv(ifi,itet).gt.w_limit) then
					vanom=v2tmp(ifi,itet)
					aver=aver+vanom
					naver=naver+1
				end if
				v2tmp(ifi,itet)=vanom
				if(kod_apriori.eq.1) then
					ttt=(itet-1)*dtmap+tmap1+tet0
					fff=(ifi-1)*dfmap+fmap1+fi0

                                    if(key_ft1_xy2.eq.1) then
                                        call SFDEC(fff,ttt,0.,xxx,yyy,Z,fi0,tet0)
                                    else
                                        xxx=fff
                                        yyy=ttt
                                    end if

					dv_aprio = vert_anom(xxx,yyy,zzz,ips)
					v2tmp(ifi,itet)=v2tmp(ifi,itet)+dv_aprio
				end if
				!if(itet.eq.101) write(*,*)vanom,vvv(ifi,itet)
			end do
		end do
		!pause
		if(naver.gt.0) then
			aver=aver/naver
		end if
		if(kod_av_bias.eq.1) v2tmp=v2tmp-aver


                open(14,file='../../../TMP_files/hor/dv'//ps//it//lv//'.xyz')
                do ifi=1,nfmap
                    fff=(ifi-1)*dfmap+fmap1+fi0
                    do itet=1,ntmap
                        ttt=(itet-1)*dtmap+tmap1+tet0
                        write(14,*)fff,ttt,v2tmp(ifi,itet)
                    end do
                end do
                close(14)


		open(14,file='../../../TMP_files/hor/dv'//ps//it//lv//'.grd')
		write(14,'(a4)')dsaa
		write(14,*)nfmap,ntmap
		write(14,*)fmap1+fi0,fmap2+fi0
		write(14,*)tmap1+tet0,tmap2+tet0
		write(14,*)-999,999
		do itet=1,ntmap
			write(14,*)(v2tmp(ifi,itet),ifi=1,nfmap)
		end do
		close(14)
		if(key_preview.eq.0) cycle


		!*********************************************************
		!*********************************************************
		!*********************************************************
		!*********************************************************
		!*********************************************************
		!*********************************************************

		open(14,file='config.txt')
		write(14,*)npix_x,npix_y
		write(14,*)'_______ Size of the picture in pixels (nx,ny)'
		write(14,*)fmap1+fi0,fmap2+fi0
		write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
		write(14,*)tmap1+tet0,tmap2+tet0
		write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
		write(14,*)tick_x,tick_y
		write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
		write(14,51)ar,md,it,ps,lv
		51 format('..\..\..\PICS\',a8,'\',a8,'\IT',a1,'\hor_dv',a1,a2,'.png')
		write(14,*)'_______ Path of the output picture'
		izzz=zzz
		if(ips.eq.1) then
			write(14,3331)	izzz
		3331 format(' P anomalies, depth=',i3,' km')
		else
			write(14,3332)	izzz
		3332 format(' S anomalies, depth=',i3,' km')
		end if
		write(14,*)'_______ Title of the plot on the upper axe'
		write(14,*)	1
		write(14,*)'_______ Number of layers'

		59 format('********************************************')

		write(14,59)
		write(14,*)	1
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,52)ps,it,lv
		52 format('..\..\..\TMP_files\hor\dv',a1,a1,a2,'.grd')
		write(14,*)'_______ Location of the GRD file'
		write(14,53)scale_line
		53 format(a20)
		write(14,*)'_______ Scale for visualization'
		write(14,*)	amp_min,amp_max
		write(14,*)'_______ scale diapason'


		open(1,file='../../../DATA/'//ar//'/map/polit_bound.bln',status='old',err=491)
		close(1)
		write(14,59)
		write(14,*)	2
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,54)ar
		54 format('..\..\..\DATA\',a8,'\map\polit_bound.bln')
		write(14,*)'_______ Location of the BLN file'
		write(14,*)	3
		write(14,*)'_______ Thickness of line in pixels'
		write(14,*)	100,0,0
		write(14,*)'_______ RGB color'
		491 continue

		open(1,file='../../../DATA/'//ar//'/map/coastal_line.bln',status='old',err=492)
		close(1)
		write(14,59)
		write(14,*)	2
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,55)ar
		55 format('..\..\..\DATA\',a8,'\map\coastal_line.bln')
		write(14,*)'_______ Location of the BLN file'
		write(14,*)	3
		write(14,*)'_______ Thickness of line in pixels'
		write(14,*)	0,0,0
		write(14,*)'_______ RGB color'
		492 continue


		write(14,59)
		write(14,*)	3
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,56)lv
		56 format('..\..\..\TMP_files\hor\ztr',a2,'.dat')
		write(14,*)'_______ Location of the DAT file'
		write(14,*)	1
		write(14,*)'_______ Symbol (1: circle, 2: square)'
		write(14,*)	5
		write(14,*)'_______ Size of dots in pixels'
		write(14,*)	250,0,0
		write(14,*)'_______ RGB color'

		close(14)


		i=system('layers.exe')


	end do

	v1tmp=-999
	do ifi=1,nfmap
		do itet=1,ntmap
			if (abs(vabs(1,ifi,itet)).gt.900..or.abs(vabs(2,ifi,itet)).gt.900.) cycle
			vpvs=vabs(1,ifi,itet)/vabs(2,ifi,itet)
			v1tmp(ifi,itet)=vpvs
		end do
	end do


	open(14,file='../../../TMP_files/hor/vpvs'//it//lv//'.grd')
	write(14,'(a4)')dsaa
	write(14,*)nfmap,ntmap
	write(14,*)fmap1+fi0,fmap2+fi0
	write(14,*)tmap1+tet0,tmap2+tet0
	write(14,*)-999,999
	do itet=1,ntmap
		write(14,*)(v1tmp(ifi,itet),ifi=1,nfmap)
	end do
	close(14)

	open(14,file='config.txt')
		write(14,*)npix_x,npix_y
		write(14,*)'_______ Size of the picture in pixels (nx,ny)'
		write(14,*)fmap1+fi0,fmap2+fi0
		write(14,*)'_______ Physical coordinates along X (xmin,xmax)'
		write(14,*)tmap1+tet0,tmap2+tet0
		write(14,*)'_______ Physical coordinates along Y (ymin,ymax)'
		write(14,*)tick_x,tick_y
		write(14,*)'_______ Spacing of ticks on axes (dx,dy)'
		write(14,41)ar,md,it,lv
		41 format('..\..\..\PICS\',a8,'\',a8,'\IT',a1,'\hor_vpvs',a2,'.png')
		write(14,*)'_______ Path of the output picture'
		izzz=zzz
		write(14,3333)	izzz
		3333 format(' Vp/Vs ratio, depth=',i3,' km')
		write(14,*)'_______ Title of the plot on the upper axe'
		write(14,*)	1
		write(14,*)'_______ Number of layers'

		write(14,59)
		write(14,*)	1
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,58)it,lv
		58 format('..\..\..\TMP_files\hor\vpvs',a1,a2,'.grd')
		write(14,*)'_______ Location of the GRD file'
		write(14,50)scale_line2
		50 format(a20)
		write(14,*)'_______ Scale for visualization'
		write(14,*)	vpvs_min,vpvs_max
		write(14,*)'_______ scale diapason'


		open(1,file='../../../DATA/'//ar//'/map/polit_bound.bln',status='old',err=495)
		close(1)
		write(14,59)
		write(14,*)	2
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,67)ar
		67 format('..\..\..\DATA\',a8,'\map\polit_bound.bln')
		write(14,*)'_______ Location of the BLN file'
		write(14,*)	3
		write(14,*)'_______ Thickness of line in pixels'
		write(14,*)	100,0,0
		write(14,*)'_______ RGB color'
		495 continue

		open(1,file='../../../DATA/'//ar//'/map/coastal_line.bln',status='old',err=496)
		close(1)
		write(14,59)
		write(14,*)	2
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,68)ar
		68 format('..\..\..\DATA\',a8,'\map\coastal_line.bln')
		write(14,*)'_______ Location of the BLN file'
		write(14,*)	3
		write(14,*)'_______ Thickness of line in pixels'
		write(14,*)	0,0,0
		write(14,*)'_______ RGB color'
		496 continue


		write(14,59)
		write(14,*)	3
		write(14,*)'_______ Key of the layer (1: contour, 2: line, 3:dots)'
		write(14,66)lv
		66 format('..\..\..\TMP_files\hor\ztr',a2,'.dat')
		write(14,*)'_______ Location of the DAT file'
		write(14,*)	1
		write(14,*)'_______ Symbol (1: circle, 2: square)'
		write(14,*)	5
		write(14,*)'_______ Size of dots in pixels'
		write(14,*)	250,0,0
		write(14,*)'_______ RGB color'

		close(14)

		i=system('layers.exe')


end do

open(11,file='../../../DATA/'//ar//'/'//md//'/data/vaver'//it//'.dat')
do ilev=1,nlev
    write(11,*)hlev(ilev),vaver(1,ilev),vaver(2,ilev)
end do
close(11)


stop
end