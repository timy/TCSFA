#define GRID_NX 300
#define GRID_NZ 600
#define GRID_LOWER_X  -2d0
#define GRID_UPPER_X   2d0 
#define GRID_LOWER_Z  -4d0
#define GRID_UPPER_Z   4d0 
#define N_LINE 50766569
#define N_TRAJ_TYPE   4

#define FID_RAW 101 
#define DATA_FILE_RAW '../../dat/data.dat'
#define FID_TRAJ 110
#define DATA_FILE_SELECT '../data/select'

#define N_SELECT 1

#define SELECT_PX   0.5d0
#define SELECT_PZ   1.5875d0

program main

    implicit none
    integer, parameter:: n_pos = N_LINE
    integer:: i_pos, i_px, i_pz, ierr, iz, ix, n_pass_x, n_pass_z
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_M_re, data_M_im
    
    
    integer, parameter:: nx = GRID_NX
    integer, parameter:: nz = GRID_NZ
    double precision, parameter:: px_upper = GRID_UPPER_X
    double precision, parameter:: px_lower = GRID_LOWER_X
    double precision, parameter:: pz_upper = GRID_UPPER_Z
    double precision, parameter:: pz_lower = GRID_LOWER_Z
    double precision, parameter:: step_length_x = ( px_upper - px_lower ) / nx;
    double precision, parameter:: step_length_z = ( pz_upper - pz_lower ) / nz;
    double precision:: grid_px(nx), grid_pz(nz), grid_w(nx,nz)
    integer:: grid_count(GRID_NX,GRID_NZ)
    integer:: ierr_read, n_err, i_type
    ! nx, nz: the number of intervals on the grid
    
    integer:: i_select,  i_select_x, i_select_z, count
    double precision:: select_px(N_SELECT), select_pz(N_SELECT)
    integer:: index_select_x(N_SELECT), index_select_z(N_SELECT)
    integer:: traj_type_count(nx,nz,N_TRAJ_TYPE)
    character(len=200):: text

    forall(ix=1:nx) grid_px(ix) = px_lower + (ix-0.5d0)*step_length_x;
    forall(iz=1:nz) grid_pz(iz) = pz_lower + (iz-0.5d0)*step_length_z;
    forall(ix=1:nx,iz=1:nz) grid_count(ix,iz) = 0;
    forall(i_px=1:nx,i_pz=1:nz,i_type=1:N_TRAJ_TYPE) &
          traj_type_count(i_px,i_pz,i_type) = 0;
    n_err = 0;

    open(FID_RAW, file=DATA_FILE_RAW, status='OLD');

    select_px = (/SELECT_PX/);
    select_pz = (/SELECT_PZ/);

    forall(i_select=1:N_SELECT) index_select_x(i_select) = &
          ceiling( nx * ( select_px(i_select) - px_lower ) / ( px_upper - px_lower ) );
    forall(i_select=1:N_SELECT) index_select_z(i_select) = &
          ceiling( nz * ( select_pz(i_select) - pz_lower ) / ( pz_upper - pz_lower ) );

    count = 0;

    do i_select = 1, N_SELECT
        write(*, '(a, 2(e15.8, 2x))'), &
              'the program will select trajectories for: ', &
              grid_px( index_select_x(i_select) ), &
              grid_pz( index_select_z(i_select) )
        write(text, '(a,i2,a)'), DATA_FILE_SELECT, i_select, '.dat';
        print*, 'the selected trajectories will be written to file: ', trim(text);
        open( FID_TRAJ+i_select, file=trim(text) );
    end do

    print*, 'start to loop the raw data ...'

    do i_pos = 1, n_pos
        
        read(FID_RAW, '(10(e15.8,1x),3(i2,1x))', err=101, iostat=ierr_read),  &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, &
              data_M_re, data_M_im, n_pass_x, n_pass_z, ierr;

        i_px = ceiling( ( data_px_inf - px_lower ) / step_length_x );
        i_pz = ceiling( ( data_pz_inf - pz_lower ) / step_length_z );

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
            i_type = 0; ! unknown
        end if


        do i_select = 1, N_SELECT
            if(i_px == index_select_x(i_select) .and. i_pz == index_select_z(i_select) ) then
                if( data_M_re**2 + data_M_im**2 > 1d-2 ) then
                    write(FID_TRAJ+i_select, '(10(e15.8,1x),4(i2,1x))'),  &
                          data_px_0, data_pz_0, data_ts_re, data_ts_im, &
                          data_x_0, data_z_0, data_px_inf, data_pz_inf, &
                          data_M_re, data_M_im, n_pass_x, n_pass_z, ierr, i_type;
                end if
            end if
        end do

        if( mod(i_pos, 2000000) == 0 ) then
            write(*, '(a,f6.2,a)'), 'progress towards completion: ', i_pos * 100d0 / N_LINE, '%' ;
        end if

        ! error
        if(ierr > 0) then
            n_err = n_err + 1;
            cycle;
        end if
        ! boundary
        if(i_px < 1 .or. i_px > nx) cycle;
        if(i_pz < 1 .or. i_pz > nz) cycle;

        grid_count(i_px,i_pz) = grid_count(i_px,i_pz) + 1;
        traj_type_count(i_px,i_pz,i_type) = traj_type_count(i_px,i_pz,i_type) + 1;
      
        cycle;

101     print*, 'error during reading line', i_pos, ':', ierr_read;
        
    end do
    close(FID_RAW);

    do i_select = 1, N_SELECT
        close( FID_TRAJ+i_select );
    end do

    print*, 'number of valid data:', N_LINE - n_err;
    print*, 'number of suspecting data:', n_err;
    print*, 'number of data in use:', sum(grid_count);
    print*, 'done!'


end program main
