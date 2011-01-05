#define PI    3.14159265
#define GRID_NX 1000
#define GRID_NZ 2000
#define GRID_N_ANGLE  10
#define GRID_N_ENERGY 250
#define GRID_LOWER_X  -2d0
#define GRID_UPPER_X  2d0 
#define GRID_LOWER_Z  -4d0
#define GRID_UPPER_Z  4d0
#define GRID_LOWER_ANGLE   0d0
#define GRID_UPPER_ANGLE   PI
#define GRID_LOWER_ENERGY  0d0
#define GRID_UPPER_ENERGY  0.5d0
#define N_LINE 63706888

#define FID_RAW 101
#define DATA_FILE_RAW '../../dat/data.dat'
#define FID_SPEC 102
#define DATA_FILE_SPEC '../data/spectra.dat'

#define N_TRAJ_TYPE 5

program main

    implicit none
    integer, parameter:: n_pos = N_LINE
    integer:: i_pos, i_px, i_pz, ierr_read, i_type, i_angle, i_energy 
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_M_re, data_M_im
    integer::  n_pass_x, n_pass_z, ierr
    
    integer, parameter:: nx = GRID_NX
    integer, parameter:: nz = GRID_NZ
    integer, parameter:: n_angle = GRID_N_ANGLE
    integer, parameter:: n_energy = GRID_N_ENERGY
    double precision, parameter:: px_upper = GRID_UPPER_X
    double precision, parameter:: px_lower = GRID_LOWER_X
    double precision, parameter:: pz_upper = GRID_UPPER_Z
    double precision, parameter:: pz_lower = GRID_LOWER_Z
    double precision, parameter:: energy_upper = GRID_UPPER_ENERGY
    double precision, parameter:: energy_lower = GRID_LOWER_ENERGY
    double precision, parameter:: angle_upper = GRID_UPPER_ANGLE
    double precision, parameter:: angle_lower = GRID_LOWER_ANGLE
    double precision, parameter:: d_px = ( px_upper - px_lower ) / nx
    double precision, parameter:: d_pz = ( pz_upper - pz_lower ) / nz
    double precision, parameter:: d_angle = (angle_upper - angle_lower) / n_angle
    double precision, parameter:: d_energy = (energy_upper - energy_lower) / n_energy

    double precision:: grid_px(nx), grid_pz(nz), grid_energy(n_energy), grid_angle(n_angle)   
    double complex:: qtm_Mp(nx,nz,N_TRAJ_TYPE)
    double precision:: cls_Mp(nx,nz,N_TRAJ_TYPE)
    character(len=200):: text
    integer:: grid_count(nx,nz), traj_type_count(nx,nz,N_TRAJ_TYPE)
    integer:: count_bound, count_error
    double precision:: angle, energy

    double precision:: w_qtm_all(n_angle,n_energy), w_qtm_1(n_angle,n_energy), w_qtm_2(n_angle,n_energy), &
          w_qtm_3(n_angle,n_energy), w_qtm_4(n_angle,n_energy)
    double precision:: w_cls_all(n_angle,n_energy), w_cls_1(n_angle,n_energy), w_cls_2(n_angle,n_energy), &
          w_cls_3(n_angle,n_energy), w_cls_4(n_angle,n_energy)


    ! initialization
    forall(i_px=1:nx) grid_px(i_px) = px_lower + (i_px-0.5d0)*d_px;
    forall(i_pz=1:nz) grid_pz(i_pz) = pz_lower + (i_pz-0.5d0)*d_pz;
    forall(i_energy=1:n_energy) grid_energy(i_energy) = energy_lower + (i_energy-0.5d0)*d_energy;
    forall(i_angle=1:n_angle) grid_angle(i_angle) = angle_lower + (i_angle-0.5d0)*d_angle;
    forall(i_px=1:nx,i_pz=1:nz) grid_count(i_px,i_pz) = 0;

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
            write(*, '(a,f6.2,a)'), 'progress towards completion: ', i_pos * 100d0 / N_LINE, '%' ;
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
        if( data_pz_inf * data_z_0 > 0 .and. data_px_inf * data_px_0 >= 0 ) then
            i_type = 1;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 > 0 ) then
            i_type = 2;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 < 0 ) then
            i_type = 3;
        elseif( data_pz_inf * data_z_0 > 0 .and. data_px_inf * data_px_0 < 0 ) then
            i_type = 4;
        elseif( data_pz_inf * data_z_0 < 0 .and. data_px_inf * data_px_0 == 0 ) then
            i_type = 5; ! HHG
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

    do i_px = 1, nx
        do i_pz = 1, nz
            
            ! identify the angle index
            angle = datan( grid_px(i_px) / grid_pz(i_pz) );

            ! map from [-pi/2, pi/2] to [0,2*pi]
            if( grid_pz(i_pz) < 0 .and. grid_px(i_px) > 0 ) then
                angle = angle + PI; ! 2nd quadrant
            end if
            if( grid_pz(i_pz) < 0 .and. grid_px(i_px) < 0 ) then
                angle = angle + PI; ! 3rd quadrant
            end if
            if( grid_pz(i_pz) > 0 .and. grid_px(i_px) < 0 ) then
                angle = angle + 2d0*PI; ! 4th quadrant
            end if

            i_angle = ceiling( (angle - angle_lower) / d_angle );
            ! identify the energy index
            energy = ( grid_px(i_px)**2 + grid_pz(i_pz)**2 ) / 2d0;
            i_energy = ceiling( ( energy - energy_lower ) / d_energy );

            ! boundary
            if(i_angle < 1 .or. i_angle > n_angle) cycle;
            if(i_energy < 1 .or. i_energy > n_energy) cycle;

            ! calculate probability in each energy bin with different amplitudes
            w_qtm_all(i_angle,i_energy) = w_qtm_all(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cdabs(sum(qtm_Mp(i_px,i_pz,1:4)))**2;
            w_qtm_1(i_angle,i_energy) = w_qtm_1(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cdabs(qtm_Mp(i_px,i_pz,1))**2;
            w_qtm_2(i_angle,i_energy) = w_qtm_2(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cdabs(qtm_Mp(i_px,i_pz,2))**2;
            w_qtm_3(i_angle,i_energy) = w_qtm_3(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cdabs(qtm_Mp(i_px,i_pz,3))**2;
            w_qtm_4(i_angle,i_energy) = w_qtm_4(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cdabs(qtm_Mp(i_px,i_pz,4))**2;
            w_cls_all(i_angle,i_energy) = w_cls_all(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * sum(cls_Mp(i_px,i_pz,1:4))**2;
            w_cls_1(i_angle,i_energy) = w_cls_1(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cls_Mp(i_px,i_pz,1)**2;
            w_cls_2(i_angle,i_energy) = w_cls_2(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cls_Mp(i_px,i_pz,2)**2;
            w_cls_3(i_angle,i_energy) = w_cls_3(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cls_Mp(i_px,i_pz,3)**2;
            w_cls_4(i_angle,i_energy) = w_cls_4(i_angle,i_energy) + &
                  dabs(grid_px(i_px)) * cls_Mp(i_px,i_pz,4)**2;
        end do
    end do


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

    print*, 'now writing the total spectra...';
    open(FID_SPEC, file='../data/energy_spec_qtm_all.dat');
    open(FID_SPEC+1, file='../data/energy_spec_qtm_1.dat');
    open(FID_SPEC+2, file='../data/energy_spec_qtm_2.dat');
    open(FID_SPEC+3, file='../data/energy_spec_qtm_3.dat');
    open(FID_SPEC+4, file='../data/energy_spec_qtm_4.dat');
    open(FID_SPEC+5, file='../data/energy_spec_cls_all.dat');
    open(FID_SPEC+6, file='../data/energy_spec_cls_1.dat');
    open(FID_SPEC+7, file='../data/energy_spec_cls_2.dat');
    open(FID_SPEC+8, file='../data/energy_spec_cls_3.dat');
    open(FID_SPEC+9, file='../data/energy_spec_cls_4.dat');

    do i_angle = 1, n_angle
        do i_energy = 1, n_energy
            
            if ( dabs(w_qtm_all(i_angle,i_energy)) < 1d-99 ) w_qtm_all(i_angle,i_energy) = 1d-99;
            write(FID_SPEC, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_qtm_all(i_angle,i_energy);
            if ( dabs(w_qtm_1(i_angle,i_energy)) < 1d-99 ) w_qtm_1(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+1, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_qtm_1(i_angle,i_energy);
            if ( dabs(w_qtm_2(i_angle,i_energy)) < 1d-99 ) w_qtm_2(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+2, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_qtm_2(i_angle,i_energy);
            if ( dabs(w_qtm_3(i_angle,i_energy)) < 1d-99 ) w_qtm_3(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+3, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_qtm_3(i_angle,i_energy);
            if ( dabs(w_qtm_4(i_angle,i_energy)) < 1d-99 ) w_qtm_4(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+4, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_qtm_4(i_angle,i_energy);

            if ( dabs(w_cls_all(i_angle,i_energy)) < 1d-99 ) w_cls_all(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+5, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_cls_all(i_angle,i_energy);
            if ( dabs(w_cls_1(i_angle,i_energy)) < 1d-99 ) w_cls_1(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+6, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_cls_1(i_angle,i_energy);
            if ( dabs(w_cls_2(i_angle,i_energy)) < 1d-99 ) w_cls_2(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+7, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_cls_2(i_angle,i_energy);
            if ( dabs(w_cls_3(i_angle,i_energy)) < 1d-99 ) w_cls_3(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+8, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_cls_3(i_angle,i_energy);
            if ( dabs(w_cls_4(i_angle,i_energy)) < 1d-99 ) w_cls_4(i_angle,i_energy) = 1d-99;
            write(FID_SPEC+9, '(3(e15.8,2x))'), grid_energy(i_energy), grid_angle(i_angle), &
                  w_cls_4(i_angle,i_energy);
                
        end do
    end do

    close(FID_SPEC);
    close(FID_SPEC+1);
    close(FID_SPEC+2);
    close(FID_SPEC+3);
    close(FID_SPEC+4);
    close(FID_SPEC+5);
    close(FID_SPEC+6);
    close(FID_SPEC+7);
    close(FID_SPEC+8);
    close(FID_SPEC+9);

    print*, 'done!'


end program main
