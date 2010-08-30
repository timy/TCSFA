#define GRID_NX 100
#define GRID_NZ 800
#define GRID_LOWER_X  0d0
#define GRID_UPPER_X  0.5d0 
#define GRID_LOWER_Z  -2d0
#define GRID_UPPER_Z  2d0 
#define N_LINE  37123188
!#define N_LINE 514

#define FID_RAW 101
#define DATA_FILE_RAW '../../dat/data.dat'
!#define DATA_FILE_RAW 'dat/select.dat'
#define FID_SPEC 102
#define DATA_FILE_SPEC '../data/spectra.dat'

#define N_TRAJ_TYPE 4
!#define N_TUNNEL_TIME 8
#define N_TUNNEL_TIME 12

! criteria to classify which cycle the traj belongs to: 1: fix interval; 2: antomative group 
#define TIME_INTERVAL_CHOOSE 2

! criteria to choose best trajectory within a family: 1: p_inf; 2: max(w)
#define CRITERA_BEST_TRAJ 1

! coordinate system: 1: cartesian; 2: cylindrical
#define SPECTRA_COORD 1

#if TIME_INTERVAL_CHOOSE == 1
#define TIME_BOUNDARY 138.9d0, 221.8d0, 274.7d0, 346.2d0, 415.8d0, 478.4d0, 549.9d0, 608.8d0, 687.4d0
#elif TIME_INTERVAL_CHOOSE == 2
#define PULSE_CYCLE 275.578d0
#endif

program main

    implicit none
    integer, parameter:: n_pos = N_LINE
    integer:: i_pos, i_px, i_pz, i_type, i_temp, ierr
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_M_re, data_M_im
    
    
    integer, parameter:: nx = GRID_NX
    integer, parameter:: nz = GRID_NZ
    double precision, parameter:: px_upper = GRID_UPPER_X
    double precision, parameter:: px_lower = GRID_LOWER_X
    double precision, parameter:: pz_upper = GRID_UPPER_Z
    double precision, parameter:: pz_lower = GRID_LOWER_Z
    double precision, parameter:: d_px = ( px_upper - px_lower ) / nx;
    double precision, parameter:: d_pz = ( pz_upper - pz_lower ) / nz;
    double precision:: grid_px(nx), grid_pz(nz), grid_w(nx,nz)
    integer:: grid_count(nx,nz)
    double complex:: grid_Mp(nx,nz)    
    integer:: ierr_read, n_pass_x, n_pass_z

    integer:: traj_type_count(nx,nz,N_TRAJ_TYPE)

    double complex:: qtm_Mp(nx,nz,N_TRAJ_TYPE)
    double precision:: cls_Mp(nx,nz,N_TRAJ_TYPE)
    double precision:: delta_px, delta_pz, r
    double precision:: dt_a, dt_b
    character(len=200):: text
    integer:: count_bound, count_error

    ! initialization
    forall(i_px=1:nx) grid_px(i_px) = px_lower + (i_px-0.5d0)*d_px;
    forall(i_pz=1:nz) grid_pz(i_pz) = pz_lower + (i_pz-0.5d0)*d_pz;
    forall(i_px=1:nx,i_pz=1:nz) grid_count(i_px,i_pz) = 0;
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px, i_pz) = (0d0, 0d0);

    forall(i_px=1:nx,i_pz=1:nz,i_type=1:N_TRAJ_TYPE) &
          traj_type_count(i_px,i_pz,i_type) = 0;

    forall(i_px=1:nx,i_pz=1:nz,i_type=1:N_TRAJ_TYPE) &
          qtm_Mp(i_px,i_pz,i_type) = (0d0,0d0);
    forall(i_px=1:nx,i_pz=1:nz,i_type=1:N_TRAJ_TYPE) &
          cls_Mp(i_px,i_pz,i_type) = 0d0;


    ! open the raw data file
    open(FID_RAW, file=DATA_FILE_RAW, status='OLD');

    print*, 'start to loop the raw data ...'

    ! start looping the file


    count_error = 0;
    count_bound = 0;
    do i_pos = 1, n_pos
        
        read(FID_RAW, '(10(e15.8,1x),3(i2,1x))', err=101, iostat=ierr_read),  &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, &
              data_M_re, data_M_im, n_pass_x, n_pass_z, ierr;

        i_px = ceiling( ( data_px_inf - px_lower ) / d_px );
        i_pz = ceiling( ( data_pz_inf - pz_lower ) / d_pz );

        if( mod(i_pos, 2000000) == 0 ) then
            write(*, '(a,f6.2,a)'), 'progress toward completetion: ', i_pos * 100d0 / N_LINE, '%' ;
        end if

        ! check whether the data is reasonable for further process
        if(ierr > 0) then
            if( ierr == 2 ) then
                count_bound = count_bound + 1;
            else
                count_error = count_error + 1;
            end if
            cycle;
        end if
        ! boundary
        if(i_px < 1 .or. i_px > nx) cycle;
        if(i_pz < 1 .or. i_pz > nz) cycle;



        grid_count(i_px, i_pz) = grid_count(i_px, i_pz) + 1;

        i_type = 0;
        ! decide which type this trajectory belongs to
        if( data_pz_inf * data_z_0 > 0 .and. data_px_inf * data_px_0 > 0 ) then
            i_type = 1;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 > 0 ) then
            i_type = 2;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 < 0 ) then
            i_type = 3;
        elseif( data_pz_inf * data_z_0 > 0 .and. data_px_inf * data_px_0 < 0 ) then
            i_type = 4;
        else
            print*, 'unknown type!';
            cycle;
        end if
        
        ! now we have known what type, and from which cycle the traj is
        traj_type_count(i_px,i_pz,i_type) = traj_type_count(i_px,i_pz,i_type) + 1;
        qtm_Mp(i_px,i_pz,i_type) = qtm_Mp(i_px,i_pz,i_type) +  dcmplx( data_M_re, data_M_im );
        cls_Mp(i_px,i_pz,i_type) = cls_Mp(i_px,i_pz,i_type) + dsqrt( data_M_re**2d0 + data_M_im**2d0 );

        cycle;

101     print*, 'error during reading line', i_pos, ':', ierr_read;
        
    end do
   
    close(FID_RAW);


    ! ==============================================================================
    ! finish the loop of the huge file, output some information

    print*, 'number of valid data in specified region:', sum(grid_count);
    print*, 'number of bound electrons:', count_bound;
    print*, 'number of divergent trajectories:', count_error;
    print*, 'number of data in use for each type:'
    do i_type = 1, N_TRAJ_TYPE
        print*, 'type', i_type, sum(traj_type_count(:,:,i_type));
    end do


    ! ==============================================================================
    ! ------------------------------------------------------------
    ! the quantum spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,:));

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_all.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);



    ! ------------------------------------------------------------
    ! the quantum spectra for each type of traj
    do i_type = 1, N_TRAJ_TYPE

        ! transition amplitude
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  qtm_Mp(i_px,i_pz,i_type);
    
        ! ionization probability
#if SPECTRA_COORD == 1
    
        print*, 'now calculate the spectra in cartesian coordinate ...'
        forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
              cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

        print*, 'now calculate the spectra in cylindrical coordinate ...'
        forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
              cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

        write(text, '(a,i2,a)'), '../data/spec_qtm_', i_type, '.dat';
        print*, 'now writing quantum spec for each type of traj ...', text;
        open(FID_SPEC, file=trim(text));
        do i_px = 1, nx
            do i_pz = 1, nz
                
                if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
                write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                      grid_w(i_px,i_pz);
                
            end do
        end do
        close(FID_SPEC);

    end do



    ! ------------------------------------------------------------
    ! the classical spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(cls_Mp(i_px,i_pz,:));

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_cls_all.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);


    ! ------------------------------------------------------------
    ! the classical spectra for each type of traj
    do i_type = 1, N_TRAJ_TYPE

        ! transition amplitude
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  cls_Mp(i_px,i_pz,i_type);
    
        ! ionization probability
#if SPECTRA_COORD == 1
    
        print*, 'now calculate the spectra in cartesian coordinate ...'
        forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
              cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

        print*, 'now calculate the spectra in cylindrical coordinate ...'
        forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
              cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

        write(text, '(a,i2,a)'), '../data/spec_cls_', i_type, '.dat';
        print*, 'now writing classical spec for each type of traj ...', text;
        open(FID_SPEC, file=trim(text));
        do i_px = 1, nx
            do i_pz = 1, nz
                
                if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
                write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                      grid_w(i_px,i_pz);
                
            end do
        end do
        close(FID_SPEC);

    end do



    ! ------------------------------------------------------------
    ! number distribution

    write(text, '(a)'), '../data/spec_number.dat';
    print*, 'now writing the spectra of number distribution ...', text
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz

            write(FID_SPEC, '(2(e15.8,2x), i8)'), grid_pz(i_pz), grid_px(i_px), &
                  grid_count(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC)    





    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type2
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,1:2));

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_1+2.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type2 + type3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,1:3));

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_123.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    ! ------------------------------------------------------------
    ! the quantum spectra for type3 + type4
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,3:4));

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_3+4.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  qtm_Mp(i_px,i_pz,1) + qtm_Mp(i_px,i_pz,3);

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_1+3.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type4
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = qtm_Mp(i_px,i_pz,1) + qtm_Mp(i_px,i_pz,4);

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_1+4.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    ! ------------------------------------------------------------
    ! the quantum spectra for type2 + type3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,2:3));

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_2+3.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    ! ------------------------------------------------------------
    ! the quantum spectra for type2 + type4
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  qtm_Mp(i_px,i_pz,2) + qtm_Mp(i_px,i_pz,4);

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2;

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'
    forall(i_px=1:nx,i_pz=1:nz) grid_w(i_px,i_pz) = &
          cdabs( grid_Mp(i_px,i_pz) )**2 * dabs( grid_px(i_px) );

#endif

    write(text, '(a)'), '../data/spec_qtm_2+4.dat';
    print*, 'now writing the total spectra...', text;
    open(FID_SPEC, file=trim(text));
    do i_px = 1, nx
        do i_pz = 1, nz
            
            if ( dabs(grid_w(i_px,i_pz)) < 1d-99 ) grid_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_pz(i_pz), grid_px(i_px), &
                  grid_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    print*, 'done!'


end program main
