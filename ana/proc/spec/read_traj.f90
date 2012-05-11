subroutine read_traj( nx, nz, grid_lower_x, d_px, &
      grid_lower_z, d_pz, n_type, rank, n_traj, &
      traj_count, err_count, qtm_Mp, cls_Mp, L )
    implicit none
    integer, parameter:: fid_traj = 102
    integer, intent(in):: nx, nz, n_type, rank, n_traj
    double precision, intent(in):: grid_lower_x, grid_lower_z, d_px, d_pz
    integer:: traj_count( nx, nz, n_type ), err_count( 4 )
    double precision:: cls_Mp( nx, nz, n_type )
    double complex:: qtm_Mp( nx, nz, n_type )
    double precision:: L( nx, nz, n_type )
    character(len=64):: file_name
    integer:: i_pos, i_px, i_pz, i_type, ierr_read, b_filter
    integer, parameter:: b_mirrow = 1

    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L1, data_L2, &
          data_M_re, data_M_im
    integer:: data_n_pass_x, data_n_pass_z, data_ierr, data_n_step
    integer, external:: filter_condition

    write( file_name, '(a,i3,a)' ), '../../../dat/traj_', rank, '.dat'
    open( fid_traj, file=file_name )

    do i_pos = 1, n_traj
        
        read( fid_traj, *, err=101, iostat=ierr_read, end=105 ), &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L1, data_L2, &
              data_M_re, data_M_im, data_n_step, data_ierr

        if( data_ierr > 0 ) then

            ! 1: cannot find saddle point
            ! 2: cannot find convergent traj for rk4
            ! 3: final energy is less than 0, bound
            ! 4: get too close to the core

            err_count( data_ierr ) = err_count( data_ierr ) + 1
            cycle

        end if

        if( b_mirrow .eq. 1 ) then
            if( data_px_inf < 0 ) then
                data_px_0 = - data_px_0
                data_px_inf = - data_px_inf
            end if
        end if

        ! identify the index of bin the current traj belongs to
        i_px = ceiling( ( data_px_inf - grid_lower_x ) / d_px );
        i_pz = ceiling( ( data_pz_inf - grid_lower_z ) / d_pz );
        
        ! check if in the spectra-ploting region
        if(i_px < 1 .or. i_px > nx) cycle;
        if(i_pz < 1 .or. i_pz > nz) cycle;

        ! decide the type of current traj
        if( data_pz_inf * data_z_0 > 0d0 .and. &
              data_px_inf * data_px_0 >= 0d0 ) then
            i_type = 1;
        elseif( data_pz_inf * data_z_0 < 0d0 .and. &
              data_px_inf * data_px_0 > 0d0 ) then
            i_type = 2;
        elseif( data_pz_inf * data_z_0 < 0d0 .and. &
              data_px_inf * data_px_0 < 0d0 ) then
            i_type = 3;
        elseif( data_pz_inf * data_z_0 > 0d0 .and. &
              data_px_inf * data_px_0 < 0d0 ) then
            i_type = 4;
        elseif( data_pz_inf * data_z_0 < 0d0 .and. &
              data_px_inf * data_px_0 == 0d0 ) then
            i_type = 5; ! HHG, moving along the z-axis, head-to-head collision
        else
            i_type = 0;
            print*, 'unknown type!';
            cycle;
        end if

        b_filter = filter_condition( data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L1, data_L2, data_M_re, data_M_im, &
              data_n_step, i_type )

        if( b_filter .eq. 1 ) then
            ! count for statistics information
            traj_count(i_px,i_pz,i_type) = traj_count(i_px,i_pz,i_type) + 1;
        
            ! superposition in "binning grid" for quantum spectra
            qtm_Mp(i_px,i_pz,i_type) = qtm_Mp(i_px,i_pz,i_type) &
                  + dcmplx( data_M_re, data_M_im );
        
            ! incoherent summation in "binning grid" for classical spectra
            cls_Mp(i_px,i_pz,i_type) = cls_Mp(i_px,i_pz,i_type) &
                  + dsqrt( data_M_re**2d0 + data_M_im**2d0 );

            ! L
            L(i_px,i_pz,i_type) = L(i_px,i_pz,i_type) &
                  + data_L1
        end if

        cycle;

101     print*, 'error during reading line', i_pos, ':', ierr_read, rank

    end do

    close( fid_traj )
    return

105 write(*, *), 'the traj_file not yet finished. EOF = ', i_pos, rank
    close( fid_traj )
    return

end subroutine read_traj



integer function filter_condition( px0, pz0, ts_re, ts_im, x0, z0, px, pz, L1, L2, M_re, M_im, n_step, i_type )
    implicit none
    double precision, intent(in):: px0, pz0, ts_re, ts_im, x0, z0, px, pz, L1, L2
    double precision:: M_re, M_im
    integer, intent(in):: n_step, i_type
    double complex, external:: pulse_A_z
    double precision:: vz0
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double precision, parameter:: PI = 2d0 * dasin(1d0)
    double precision, parameter:: phi = 1.0d0 * PI
    double complex:: Mp

    filter_condition = 1

!    if( M_re*M_re + M_im*M_im > 0.6**2) filter_condition = 0

    if( ( i_type .eq. 3 ) .or. ( i_type .eq. 4 ) ) then
        !         Mp = dcmplx( M_re, M_im )
        Mp = Mp * 0.2d0
        !         M_re = dreal( Mp )
        !         M_im = dimag( Mp )
    end if

!!$    if( (i_type .eq. 2) .and. n_step > 0 ) then
!!$        filter_condition = 0
!!$    end if
!!$
   
!!$    vz0 = pz0 + dreal( pulse_A_z( dcmplx(ts_re, 0d0) ) )
!!$    if( ( vz0 * z0 ) .lt. 0d0 ) then
!!$        filter_condition = 0
!!$    end if

    return
end function filter_condition


#include "../../../src/include/inc_field.h"
double complex function pulse_A_z( t ) result(A_z)
    implicit none;
    double complex, intent(in):: t;
    double precision, parameter:: Ef0 = E0
    double precision, parameter:: om = OM
    double precision, parameter:: xi = XI
    
    A_z = -( ( E0 * cdsin(om*t) ) / ( om * dsqrt( 1d0 + xi**2 ) ) )

    return;
end function pulse_A_z
