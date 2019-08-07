USE DFPORT

character*8 ar,ar_all(100),md,md_all(100),line
character*1 it
integer kod_loc(100),kod_iter(100),kod_oe(100),kod_rf(100)

i=system('mkdir ..\..\..\TMP_files\tmp')
i=system('mkdir ..\..\..\TMP_files\hor')
i=system('mkdir ..\..\..\TMP_files\rays')
i=system('mkdir ..\..\..\TMP_files\1D_mod')
i=system('mkdir ..\..\..\TMP_files\vert')

open(1, file='../../../all_areas.dat')
do i=1,4
	read(1,*)
end do
do i=1,10
	read(1,'(a8,1x,a8,1x,i1)',end=7)ar_all(i),md_all(i),kod_iter(i)
end do	
7 close(1)
n_ar=i-1
write(*,*)' n_ar=',n_ar

do iar=1,n_ar
    ar=ar_all(iar)
    md=md_all(iar)
    niter=kod_iter(iar)

    kod_param=1
    key_ft1_xy2=1
    open(1,file='../../../data/'//ar//'/'//md//'/MAJOR_PARAM.DAT',status='old',err=562)
    do i=1,10000
        read(1,'(a8)',end=563)line
        if(line.eq.'GENERAL ') goto 564
    end do
    562 continue
    write(*,*)' file MAJOR_PARAM.DAT does not exist in ar=',ar,' md=',md
    stop
    563 continue
    write(*,*)' cannot find GENERAL INFORMATION in MAJOR_PARAM.DAT!!!'
    pause
    564 continue
    read(1,*)key_1real_2syn
    read(1,*)VPSX_key
    read(1,*)koe
    read(1,*)kref
    read(1,*,err=565)key_ft1_xy2
    565 close(1)


    key_table1_line2=1
    open(1,file='../../../data/'//ar//'/'//md//'/MAJOR_PARAM.DAT',status='old',err=562)
    do i=1,10000
        read(1,'(a8)',end=585)line
        if(line.eq.'1D LOCAT') goto 584
    end do
    582 continue
    write(*,*)' file MAJOR_PARAM.DAT does not exist in ar=',ar,' md=',md
    pause
    584 continue
    read(1,*)key_table1_line2
    585 close(1)

    open(11,file='../../../model.dat')
    write(11,'(a8)')ar		
    write(11,'(a8)')md		
    close(11)

    i=system('mkdir ..\..\..\DATA\'//ar//'\'//md//'\data')

    write(*,*)' DATASET:',ar,' MODEL:',md
    write(*,*)' key_ft1_xy2=',key_ft1_xy2,' kref=',kref

    !GOTO 771
    if(key_1real_2syn.eq.2) then
         write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZATION OF THE SYNTHETIC MODEL IN HORIZONTAL SECTIONS'
        i=system('..\..\4_CREATE_SYN_DATA\a_set_syn_hor\create.exe')
         write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZATION OF THE SYNTHETIC MODEL IN VERTICAL SECTIONS'
        i=system('..\..\4_CREATE_SYN_DATA\a_set_syn_ver\create.exe')

        write(*,*)' _________________________________________________________'
        write(*,*)' COMPUTE SYNTHETIC TRAVEL TIMES'
       if (key_ft1_xy2.eq.1) then
            i=system('..\..\4_CREATE_SYN_DATA\b_synth_times\rays.exe')
        else
            i=system('..\..\4_CREATE_SYN_DATA\b_xyz_synth_times\rays.exe')
        end if
    end if

    i=system('copy ..\..\..\DATA\'//ar//'\'//md//'\ref_start.dat &
    ..\..\..\DATA\'//ar//'\'//md//'\data\refmod.dat')

    if(kref.eq.1) then
! Optimization for the reference model (only in case of geographical coordinates)
        if (key_ft1_xy2.eq.2) then
            write(*,*)' Optimization of the reference model is provided only for'
            write(*,*)' the case of geographical coordinates.'
            write(*,*)' For the XYZ coordinates, this option is not ready yet. Sorry.'
            stop
        end if
               
        if (key_table1_line2.eq.1) then
            write(*,*)' _________________________________________________________'
            write(*,*)' COMPUTING THE REFERENCE TABLE WITH THE STARTING 1D MODEL'
            i=system('..\..\1_PRELIM_LOC\0_1_ref_rays\refrays.exe')
             write(*,*)' _________________________________________________________'
            write(*,*)' LOCALIZATION OF SOURCES USING THE 1D REFERENCE TABLE'
            i=system('..\..\1_PRELIM_LOC\0_2_loc_event\locate.exe')
             write(*,*)' _________________________________________________________'
            write(*,*)' PERFORM THE OPTIMISATION FOR THE 1D VELOCITY MODEL'
            i=system('..\..\1_PRELIM_LOC\START_1D\start_real.exe')
        else
            write(*,*)' _________________________________________________________'
            write(*,*)' PRELIMINARY LOCALIZATION OF SOURCES USING STRAIGHT RAYS'
            i=system('..\..\1_PRELIM_LOC\loc_straight\loc_1D.exe')
            write(*,*)' _________________________________________________________'
            write(*,*)' PERFORM THE OPTIMISATION FOR THE 1D VELOCITY MODEL'
            i=system('..\..\1_PRELIM_LOC\_START_1D_line\start_real.exe')
       end if
    end if
    if(niter.eq.0) cycle
         
    if (key_table1_line2.eq.1) then
        write(*,*)' _________________________________________________________'
        write(*,*)' COMPUTING THE REFERENCE TABLE WITH THE UPDATED 1D MODEL'
        i=system('..\..\1_PRELIM_LOC\0_1_ref_rays\refrays.exe')
         write(*,*)' _________________________________________________________'
        write(*,*)' LOCALIZATION OF SOURCES USING THE 1D REFERENCE TABLE'
       i=system('..\..\1_PRELIM_LOC\0_2_loc_event\locate.exe')
    else
        write(*,*)' _________________________________________________________'
        write(*,*)' PRELIMINARY LOCALIZATION OF SOURCES USING STRAIGHT RAYS'
        i=system('..\..\1_PRELIM_LOC\loc_straight\loc_1D.exe')
    end if

    771 continue

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
    kod_param=1
    open(1,file='../../../DATA/'//ar//'/'//md//'/MAJOR_PARAM.DAT')
    do i=1,10000
        read(1,'(a8)',end=513)line
        if(line.eq.'GRID_PAR') goto 514
    end do
    513 continue
    write(*,*)' cannot find GRID PARAMETERS in MAJOR_PARAM.DAT!!!'
    pause
    514 continue
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)kod_param
    close(1)
    !******************************************************************

    kod_param=1
    ! Execute the ITERATIONS:
    do iter=1,niter	
        write(it,'(i1)')iter
        open(11,file='../../../model.dat')
        write(11,'(a8)')ar		
        write(11,'(a8)')md		
        write(11,'(i1)')iter		
        close(11)

        write(*,*)' _________________________________________________________'
        write(*,*)' LOCATE THE SOURCES USING THE 3D RAY TRACING '
        i=system('..\..\2_INVERS_3D\1_locate\locate.exe')

        do igr=1,nornt
            open(11,file='../../../model.dat')
            write(11,'(a8)')ar		
            write(11,'(a8)')md		
            write(11,'(i1)')iter		
            write(11,'(i1)')igr	
            close(11)

            if(iter.eq.1) then
                write(*,*)' _________________________________________________________'
                write(*,*)' COMPUTE THE RAY DENSITY '
                i=system('..\..\2_INVERS_3D\2n_ray_density\plotray.exe')
                write(*,*)' _________________________________________________________'
                write(*,*)' DEFINE THE PARAMETERIZATION GRID '
                i=system('..\..\2_INVERS_3D\3n_grid\grid.exe')
                i=system('..\..\2_INVERS_3D\4_links\add_matr.exe')
                !write(*,*)' _________________________________________________________'
                !write(*,*)' VISUALIZE THE RAY PATHS AND GRID IN HORIZONTAL AND VERTICAL SECTIONS '
                !i=system('..\..\3_VISUAL\_vis_ray_path\paths.exe')
            end if

            write(*,*)' _________________________________________________________'
            write(*,*)' COMPUTE THE 1ST DERIVATIVE MATRIX (VP-VS SCHEME) '
            i=system('..\..\2_INVERS_3D\5_matrix\matr.exe')
            write(*,*)' _________________________________________________________'
            write(*,*)' PERFORM THE INVERSION (VP-VS SCHEME) '
            i=system('..\..\2_INVERS_3D\6_inversion\Invbig.exe')
        end do

        write(*,*)' _________________________________________________________'
        write(*,*)' COMPUTE THE VELOCITY FIELS IN 3D REGILAR GRID (VP-VS SCHEME) '
        i=system('..\..\2_INVERS_3D\7_3D_model\mod_3D.exe')
        write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZE THE RESULT IN HORIZONTAL SECTIONS (VP-VS SCHEME) '
        i=system('..\..\3_VISUAL\_vis_n_hor_result\visual.exe')
        write(*,*)' _________________________________________________________'
        write(*,*)' VISUALIZE THE RESULT IN VERTICAL SECTIONS (VP-VS SCHEME) '
        i=system('..\..\3_VISUAL\_vis_n_ver_result\visual.exe')
    end do ! Different iterations
    write(*,*)' _________________________________________________________'
    write(*,*)' CREATING THE REPORT ABOUT THE VARIANCE REDUCTION '
    i=system('..\..\3_VISUAL\variance_reduction\var_red.exe')
end do	! Different areas

stop
end