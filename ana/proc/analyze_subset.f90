#define FIG_NX 100
#define FIG_NZ 400
#define GRID_NX 100
#define GRID_NZ 400
#define GRID_LOWER_X  0d0
#define GRID_UPPER_X  0.5d0 
#define GRID_LOWER_Z  -1d0
#define GRID_UPPER_Z  1d0 

#define N_LINE  145122

#define FID_RAW 101
#define DATA_FILE_RAW '../data/dat_subset.dat'
#define FID_SPEC 102
#define OUTPUT_DIR '../data/subset/'
#define FID_T0_DIST 103
#define FID_L_DIST 104
#define FID_Z0_DIST 105
#define FID_SUBSET 106

#define N_TRAJ_TYPE 5

! coordinate system: 1: cartesian; 2: cylindrical
#define SPECTRA_COORD 1

! the following parameters are used for statistic information in specified area
#define N_CYCLE 3
#define OMEGA 0.0227817
#define N_BIN_PER_CYCLE 8000
! tunnel time bin
#define PI ( 2d0*asin(1d0) )
#define N_TIME_BIN (N_CYCLE*N_BIN_PER_CYCLE)
#define TIME_END (N_CYCLE * 2 * PI / OMEGA)
! angular momentum bin
#define N_L_BIN 3200
#define RANGE_L_A -20d0
#define RANGE_L_B 30d0
! tunnel exit bin
#define N_Z0_BIN 800
#define RANGE_Z0_A -30d0
#define RANGE_Z0_B 30d0

#define IF_OUTPUT_SUBSET .true.

program main

    implicit none
    integer, parameter:: nx = GRID_NX
    integer, parameter:: nz = GRID_NZ
    integer:: i_pos, i_px, i_pz, ierr_read, i_type, i_t0, i_L, i_z0
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, data_L, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_M_re, data_M_im
    integer::  data_n_pass_x, data_n_pass_z, data_ierr
    double precision, parameter:: d_px = ( GRID_UPPER_X - GRID_LOWER_X ) / GRID_NX
    double precision, parameter:: d_pz = ( GRID_UPPER_Z - GRID_LOWER_Z ) / GRID_NZ
    double precision, parameter:: dt = ( TIME_END ) / N_TIME_BIN
    double precision, parameter:: dL = ( RANGE_L_B - RANGE_L_A ) / N_L_BIN
    double precision, parameter:: d_z0 = ( RANGE_Z0_B - RANGE_Z0_A ) / N_Z0_BIN
    double complex:: qtm_Mp(nx,nz,N_TRAJ_TYPE)
    double precision:: cls_Mp(nx,nz,N_TRAJ_TYPE)
    double complex:: grid_Mp(nx,nz)
    character(len=200):: file_name
    integer:: traj_count(nx,nz,N_TRAJ_TYPE)
    integer:: count_bound, count_error, count_caught, count_subset
    double complex:: t0_weight(N_TIME_BIN,N_TRAJ_TYPE), L_weight(N_L_BIN,N_TRAJ_TYPE), &
          z0_weight(N_Z0_BIN,N_TRAJ_TYPE)
    double precision:: L_max, L_min, z0_max, z0_min
    logical, external:: filter_t0, filter_type, filter_L

    ! initialization
    forall(i_px=1:nx, i_pz=1:nz, i_type=1:N_TRAJ_TYPE) traj_count(i_px,i_pz,i_type) = 0;
    forall(i_px=1:nx, i_pz=1:nz, i_type=1:N_TRAJ_TYPE) qtm_Mp(i_px,i_pz,i_type) = (0d0, 0d0);
    forall(i_px=1:nx, i_pz=1:nz, i_type=1:N_TRAJ_TYPE) cls_Mp(i_px,i_pz,i_type) = 0d0;
    forall(i_t0=1:N_TIME_BIN, i_type=1:N_TRAJ_TYPE) t0_weight(i_t0,i_type) = (0d0, 0d0);
    forall(i_L=1:N_L_BIN, i_type=1:N_TRAJ_TYPE) L_weight(i_L,i_type) = (0d0, 0d0);
    forall(i_z0=1:N_Z0_BIN, i_type=1:N_TRAJ_TYPE) z0_weight(i_z0,i_type) = (0d0, 0d0);
    count_error = 0; count_bound = 0; count_caught = 0; count_subset = 0;
    L_max = 0d0; L_min = 0d0; z0_max = 0d0; z0_min = 0d0;

    ! open file of raw data
    open(FID_RAW, file=DATA_FILE_RAW, status='OLD');
    if( IF_OUTPUT_SUBSET ) open( FID_SUBSET, file=OUTPUT_DIR // 'dat_subset.dat')

    write(*,*), 'start to loop the raw data ...'
    
    ! start the loop reading lines one by one
    do i_pos = 1, N_LINE
        
        read(FID_RAW, '(11(e15.8,1x),3(i2,1x))', err=101, iostat=ierr_read),  &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
              data_M_re, data_M_im, data_n_pass_x, data_n_pass_z, &
              data_ierr;

        ! identify the index of bin the current traj belongs to
        i_px = ceiling( ( data_px_inf - GRID_LOWER_X ) / d_px );
        i_pz = ceiling( ( data_pz_inf - GRID_LOWER_Z ) / d_pz );
        i_t0 = ceiling( ( data_ts_re - 0d0 ) / dt );
        i_L  = ceiling( ( data_L - RANGE_L_A ) / dL );
        i_z0 = ceiling( ( data_z_0 - RANGE_Z0_A ) / d_z0);
        ! print information for current progress
        if( mod(i_pos, 2000000) == 0 ) then
            write(*, '(a,f6.2,a)'), 'progress towards completion: ', i_pos * 100d0 / N_LINE, '%' ;
        end if

        ! check the validity of the data
        if( data_ierr > 0 ) then
            if( data_ierr == 2 ) then           ! final energy is negative - does not ionize
                count_bound = count_bound + 1;
            else if ( data_ierr == 4 ) then     ! get to close to nuclei, high order effects
                count_caught = count_caught + 1;
            else                                ! simply something goes wrong
                count_error = count_error + 1;
            end if

            cycle;
        end if

        ! check if in the spectra-ploting region
        if(i_px < 1 .or. i_px > nx) cycle;
        if(i_pz < 1 .or. i_pz > nz) cycle;

        ! decide the type of current traj
        if( data_pz_inf * data_z_0 > 0 .and. data_px_inf * data_px_0 >= 0 ) then
            i_type = 1;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 > 0 ) then
            i_type = 2;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 < 0 ) then
            i_type = 3;
        elseif( data_pz_inf * data_z_0 > 0 .and. data_px_inf * data_px_0 < 0 ) then
            i_type = 4;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 == 0 ) then
            i_type = 5; ! HHG, moving along the z-axis, head-to-head collision
        else
            i_type = 0;
            print*, 'unknown type!';
            cycle;
        end if
        
        ! count for statistics information
        traj_count(i_px,i_pz,i_type) = traj_count(i_px,i_pz,i_type) + 1;

        ! filter:
        ! ------------------------------------------------------------------------------
        if( .not. ( filter_t0(data_ts_re) .and. filter_type(i_type) .and. filter_L(data_L) ) ) then        
            cycle;
        end if
        ! ------------------------------------------------------------------------------

        t0_weight(i_t0,i_type) = t0_weight(i_t0,i_type) + dcmplx( data_M_re, data_M_im );
        L_weight(i_L,i_type) = L_weight(i_L,i_type) + dcmplx( data_M_re, data_M_im );
        z0_weight(i_z0,i_type) = z0_weight(i_z0,i_type) + dcmplx( data_M_re, data_M_im );
        if(L_max < data_L) L_max = data_L;
        if(L_min > data_L) L_min = data_L;
        if(z0_max < data_z_0) z0_max = data_z_0;
        if(z0_min > data_z_0) z0_min = data_z_0;
        if( IF_OUTPUT_SUBSET ) then
            write(FID_SUBSET, '(11(e15.8,1x),3(i2,1x))'),  &
                  data_px_0, data_pz_0, data_ts_re, data_ts_im, &
                  data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
                  data_M_re, data_M_im, data_n_pass_x, data_n_pass_z, &
                  data_ierr;
        end if
        count_subset = count_subset + 1;
        ! superposition in "binning grid" for quantum spectra
        qtm_Mp(i_px,i_pz,i_type) = qtm_Mp(i_px,i_pz,i_type) +  dcmplx( data_M_re, data_M_im );
        
        ! incoherent summation in "binning grid" for classical spectra
        cls_Mp(i_px,i_pz,i_type) = cls_Mp(i_px,i_pz,i_type) + dsqrt( data_M_re**2d0 + data_M_im**2d0 );

        cycle;

101     print*, 'error during reading line', i_pos, ':', ierr_read;
        
    end do
   
    close(FID_RAW);
    if( IF_OUTPUT_SUBSET ) close( FID_SUBSET );
    ! ==============================================================================
    ! finish the loop of raw-data file, output statistical information

    write(*, '(a)'), 'In the spectra-ploting region:'
    write(*, '(a,2x,i12)'), 'count_valid = ', sum( traj_count );
    write(*, '(a,2x,i12)'), 'count_bound = ', count_bound;
    write(*, '(a,2x,i12)'), 'count_catch = ', count_caught;
    write(*, '(a,2x,i12)'), 'count_error = ', count_error;
    write(*, '(a)'), 'number of data in use for each type:'
    do i_type = 1, N_TRAJ_TYPE
        write(*, '(a,i2,a,2x,i12)'), 'count_type', i_type, ':', sum( traj_count(:,:,i_type) );
    end do
    write(*, '(2(a,2x,f15.8))'), 'L_min =', L_min, 'L_max =', L_max;
    write(*, '(2(a,2x,f15.8))'), 'z0_min =', z0_min, 'z0_max =', z0_max;
    if( IF_OUTPUT_SUBSET ) then
        write(*, '(a,2x,i12)'), 'count_subset = ', count_subset;
    end if    
    ! ==============================================================================
    ! ------------------------------------------------------------
    ! the quantum spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum( qtm_Mp(i_px,i_pz,1:4) );
    file_name = OUTPUT_DIR // 'spec_qtm_all.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the classical spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum( cls_Mp(i_px,i_pz,1:4) );
    file_name = OUTPUT_DIR // 'spec_cls_all.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for each type of traj
    do i_type = 1, N_TRAJ_TYPE

        ! transition amplitude
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  qtm_Mp(i_px,i_pz,i_type);
        write(file_name, '(a,i2,a)'), OUTPUT_DIR // 'spec_qtm_', i_type, '.dat';
        call w_from_grid_to_fig( grid_Mp, file_name );

    end do

    ! ------------------------------------------------------------
    ! the classical spectra for each type of traj
    do i_type = 1, N_TRAJ_TYPE

        ! transition amplitude
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  cls_Mp(i_px,i_pz,i_type);
        write(file_name, '(a,i2,a)'), OUTPUT_DIR // 'spec_cls_', i_type, '.dat';
        call w_from_grid_to_fig( grid_Mp, file_name );

    end do

    ! ------------------------------------------------------------
    ! number distribution
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = dcmplx( dble( sum( traj_count(i_px,i_pz,1:4) ) ), 0d0 );
    write(file_name, '(a)'), OUTPUT_DIR // 'spec_number.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type2
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,1:2));
    file_name = OUTPUT_DIR // 'spec_qtm_1+2.dat'
    call w_from_grid_to_fig( grid_Mp, file_name )

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type2 + type3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,1:3));
    file_name = OUTPUT_DIR // 'spec_qtm_123.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for type3 + type4
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,3:4));
    file_name = OUTPUT_DIR // 'spec_qtm_3+4.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  qtm_Mp(i_px,i_pz,1) + qtm_Mp(i_px,i_pz,3);
    file_name = OUTPUT_DIR // 'spec_qtm_1+3.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for type1 + type4
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = qtm_Mp(i_px,i_pz,1) + qtm_Mp(i_px,i_pz,4);
    file_name = OUTPUT_DIR // 'spec_qtm_1+4.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for type2 + type3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  sum(qtm_Mp(i_px,i_pz,2:3));
    file_name = OUTPUT_DIR // 'spec_qtm_2+3.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    ! ------------------------------------------------------------
    ! the quantum spectra for type2 + type4
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) =  qtm_Mp(i_px,i_pz,2) + qtm_Mp(i_px,i_pz,4);
    file_name = OUTPUT_DIR // 'spec_qtm_2+4.dat';
    call w_from_grid_to_fig( grid_Mp, file_name );

    call output_t0_distribution( t0_weight );
    call output_L_distribution( L_weight );
    call output_z0_distribution( z0_weight );
    print*, 'done!'


end program main



subroutine w_from_grid_to_fig( grid_Mp, file_name )
    implicit none

    double precision, parameter:: d_px = ( GRID_UPPER_X - GRID_LOWER_X ) / GRID_NX
    double precision, parameter:: d_pz = ( GRID_UPPER_Z - GRID_LOWER_Z ) / GRID_NZ
    double complex, intent(in):: grid_Mp( GRID_NX, GRID_NZ )
    double precision:: fig_px(FIG_NX), fig_pz(FIG_NZ), fig_w(FIG_NX, FIG_NZ )
    character(len=200):: file_name
    integer:: i_px, i_pz, ni_x, ni_z

    forall(i_px=1:FIG_NX) fig_px(i_px) = GRID_LOWER_X + (i_px - 0.5d0) * d_px;
    forall(i_pz=1:FIG_NZ) fig_pz(i_pz) = GRID_LOWER_Z + (i_pz - 0.5d0) * d_pz;
    forall(i_px=1:FIG_NX,i_pz=1:FIG_NZ) fig_w(i_px,i_pz) = 0d0;
    ni_x = GRID_NX / FIG_NX;
    ni_z = GRID_NZ / FIG_NZ;

#if SPECTRA_COORD == 1
    
    print*, 'now calculate the spectra in cartesian coordinate ...'

    do i_px = 1, FIG_NX
        do i_pz = 1, FIG_NZ
            fig_w(i_px,i_pz) = fig_w(i_px,i_pz) + &
                  sum( cdabs( grid_Mp(((i_px-1)*ni_x+1):(i_px*ni_x),((i_pz-1)*ni_z+1):(i_pz*ni_z)) )**2 );
        end do
    end do

#elif SPECTRA_COORD == 2

    print*, 'now calculate the spectra in cylindrical coordinate ...'

    do i_px = 1, FIG_NX
        do i_pz = 1, FIG_NZ
            fig_w(i_px,i_pz) = fig_w(i_px,i_pz) + &
                  sum( cdabs( grid_Mp(((i_px-1)*ni_x+1):(i_px*ni_x),((i_pz-1)*ni_z+1):(i_pz*ni_z)) )**2 );
        end do
        fig_w(i_px,i_pz) = fig_w(i_px,i_pz) * dabs( fig_px(i_px) );
    end do

#endif


    ! now write into file
    write(*, '(a, a)'), 'now writing to file: ', file_name;
    open(FID_SPEC, file=trim(file_name));
    do i_px = 1, FIG_NX
        do i_pz = 1, FIG_NZ
            
            if ( dabs(fig_w(i_px,i_pz)) < 1d-99 ) fig_w(i_px,i_pz) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), fig_pz(i_pz), fig_px(i_px), fig_w(i_px,i_pz);
                
        end do
    end do
    close(FID_SPEC);

    return;
end subroutine w_from_grid_to_fig


subroutine output_t0_distribution( t0_weight )
    implicit none
    double complex, intent(in):: t0_weight(N_TIME_BIN,N_TRAJ_TYPE)
    double precision:: t0
    double precision, parameter:: dt = ( TIME_END ) / N_TIME_BIN
    integer:: i_t0, i_type
    character(len=200):: file_name

    write(file_name, '(a)'), OUTPUT_DIR // 't0_dist.dat';
    open(FID_T0_DIST, file=trim(file_name));
    do i_t0 = 1, N_TIME_BIN
        t0 = (i_t0 - 0.5d0) * dt;
        write(FID_T0_DIST, '(6(e15.8,2x))'), t0, ( cdabs(t0_weight(i_t0, i_type)), i_type=1,N_TRAJ_TYPE ); 
    end do
    close(FID_T0_DIST);

    return;
end subroutine output_t0_distribution

subroutine output_L_distribution( L_weight )
    implicit none
    double complex, intent(in):: L_weight(N_L_BIN,N_TRAJ_TYPE)
    double precision:: L
    double precision, parameter:: dL = ( RANGE_L_B - RANGE_L_A ) / N_L_BIN
    integer:: i_L, i_type
    character(len=200):: file_name
    
    write(file_name, '(a)'), OUTPUT_DIR // 'L_dist.dat';
    open(FID_L_DIST, file=trim(file_name));
    do i_L = 1, N_L_BIN
        L = RANGE_L_A + (i_L - 0.5d0) * dL;
        write(FID_L_DIST, '(6(e15.8,2x))'), L, ( cdabs(L_weight(i_L, i_type)), i_type=1,N_TRAJ_TYPE );
    end do
    close(FID_L_DIST);

end subroutine output_L_distribution

subroutine output_z0_distribution( z0_weight )
    implicit none
    double complex, intent(in):: z0_weight(N_Z0_BIN,N_TRAJ_TYPE)
    double precision:: z0
    double precision, parameter:: d_z0 = ( RANGE_Z0_B - RANGE_Z0_A ) / N_Z0_BIN
    integer:: i_z0, i_type
    character(len=200):: file_name
    
    write(file_name, '(a)'), OUTPUT_DIR // 'z0_dist.dat';
    open(FID_Z0_DIST, file=trim(file_name));
    do i_z0 = 1, N_Z0_BIN
        z0 = RANGE_Z0_A + (i_z0 - 0.5d0) * d_z0;
        write(FID_Z0_DIST, '(6(e15.8,2x))'), z0, ( cdabs(z0_weight(i_z0, i_type)), i_type=1,N_TRAJ_TYPE );
    end do
    close(FID_Z0_DIST);
    
end subroutine output_z0_distribution

logical function filter_t0( t0 )
    
    implicit none
    double precision, intent(in):: t0

    if( t0 > 405d0 .and. t0 < 420d0 ) then
        filter_t0 = .true.;    
    else
        filter_t0 = .false.;
    end if

    return;

end function filter_t0

logical function filter_L( L )
    
    implicit none
    double precision, intent(in):: L

  !  if( ( L > 20d0 .and. L < 30d0 ) .or. L < 10d0) then ! (L > -4.72782d0 .and. L < 4d0) .or. 
        filter_L = .true.;
  !  else 
  !      filter_L = .false.;
  !  end if

    return;

end function filter_L

logical function filter_type( itype )
    
    implicit none
    integer, intent(in):: itype
    
    filter_type = .true.;

    return;

end function filter_type
