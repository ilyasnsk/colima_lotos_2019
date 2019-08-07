USE DFPORT
character*8 ar,ar_all(100),md,md_all(100),line
character*1 it
integer kod_loc(100),kod_iter(100),kod_oe(100),kod_rf(100)

i=system('mkdir ..\..\..\TMP_files\tmp')
i=system('mkdir ..\..\..\TMP_files\hor')
i=system('mkdir ..\..\..\TMP_files\rays')
i=system('mkdir ..\..\..\TMP_files\vert')

open(1, file='../../../all_areas.dat')
do i=1,4
	read(1,*)
end do
do i=1,10
	read(1,'(a8,1x,a8,1x,i1)',end=7)ar_all(i),md_all(i)
end do	
7 close(1)
n_ar=i-1
write(*,*)' n_ar=',n_ar

do iar=1,n_ar
    ar=ar_all(iar)
    md=md_all(iar)

    write(*,*)' ar=',ar,' md=',md


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
    write(*,*)' Number of grids=',nornt

    !******************************************************************
    kod_param=1
    open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
    do i=1,10000
        read(1,'(a8)',end=513)line
        if(line.eq.'ATTENUAT') goto 514
    end do
    513 continue
    write(*,*)' cannot find GRID PARAMETERS in MAJOR_PARAM.DAT!!!'
    pause
    514 continue
    read(1,*)npmin,nsmin
    read(1,*)iter
    close(1)
    !******************************************************************

    write(it,'(i1)')iter
    open(11,file='../../../model.dat')
    write(11,'(a8)')ar		
    write(11,'(a8)')md		
    close(11)

    write(*,*)' _________________________________________________________'
    write(*,*)' RAY TRACING AND READING ATTENUATION DATA '
    i=system('..\..\5_ATTENUATION\1_tracing\tracing.exe')

    open(1,file='../../../DATA/'//ar//'/'//md//'/data/num_ps_att.dat')
    read(1,*)np,ns
    close(1)

    do ips=1,2
        if(ips.eq.1.and.np.lt.npmin) then
            write(*,*)' Number of P-data is not sufficient: np=',np
            cycle
        end if
        
        if(ips.eq.2.and.ns.lt.nsmin) then
            write(*,*)' Number of S-data is not sufficient: ns=',nsp
            cycle
        end if
    
        do igr=1,nornt
            open(11,file='../../../model.dat')
            write(11,'(a8)')ar		
            write(11,'(a8)')md		
            write(11,'(i1)')ips		
            write(11,'(i1)')igr	
            close(11)

            write(*,*)' _________________________________________________________'
            write(*,*)' COMPUTE THE RAY DENSITY '
            i=system('..\..\5_ATTENUATION\2_ray_density\plotray.exe')
            write(*,*)' _________________________________________________________'
            write(*,*)' DEFINE THE PARAMETERIZATION GRID '
            i=system('..\..\5_ATTENUATION\3_grid\grid.exe')
            i=system('..\..\5_ATTENUATION\4_tetrad\Tetrad.exe')
            i=system('..\..\5_ATTENUATION\5_sosedi\add_matr.exe')
            write(*,*)' _________________________________________________________'
            write(*,*)' MATRIX CALCULATION '
            i=system('..\..\5_ATTENUATION\6_matr\matr.exe')
            write(*,*)' _________________________________________________________'
            write(*,*)' MATRIX CALCULATION '
            i=system('..\..\5_ATTENUATION\7_invers\Invbig.exe')
        end do

        write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZE THE RESULT IN HORIZONTAL SECTIONS  '
        i=system('..\..\5_ATTENUATION\8_visual_hor\visual.exe')
        write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZE THE RESULT IN VERTICAL SECTIONS  '
        i=system('..\..\5_ATTENUATION\8_visual_ver\visual.exe')
    end do
end do	! Different areas

stop
end