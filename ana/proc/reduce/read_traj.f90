subroutine read_traj( rank, n_traj, fid_output, selc_count )

    implicit none
    integer, intent(in):: fid_output
    integer, parameter:: fid_traj = 102
    integer, intent(in):: rank, n_traj
    integer, intent(out):: selc_count
    character(len=64):: filename
    integer:: err_count(4)
    integer:: i_pos, i_type, ierr_read
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
          data_M_re, data_M_im
    integer:: data_n_near_core, data_ierr
    integer, parameter:: b_mirrow = 1

    err_count = 0
    selc_count = 0
    ! open file of raw data
    write( filename, '(a,i3,a)' ), '../../../dat/traj_', rank, '.dat'
    open( fid_traj, file=filename )

    do i_pos = 1, n_traj

        ! read in raw data
        read( fid_traj, '(11(e16.8,1x),i4,1x,i1)', &
              err=101, iostat=ierr_read, end=105 ), &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
              data_M_re, data_M_im, data_n_near_core, data_ierr;

        ! check errors
        if( data_ierr > 0 ) then
            err_count( data_ierr ) = err_count(data_ierr) + 1
            cycle
        end if

        ! some trajectory selection criteria
        if( data_M_re*data_M_re + data_M_im*data_M_im < 1d-14 ) then
            cycle
        end if

        if( b_mirrow .eq. 1 ) then
            if( data_px_inf < 0 ) then
                data_px_inf = - data_px_inf
                data_px_0 = - data_px_0
            end if
        end if

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
        end if

        write( fid_output, '(11(e16.8,1x),i4,1x,i1)' ), &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
              data_M_re, data_M_im, data_n_near_core, i_type
        selc_count = selc_count + 1

        cycle

101     print*, 'error during reading line', i_pos, ':', ierr_read;
        
    end do

    ! close files
    close( fid_traj )
    return

105 write(*, *), 'the traj_file not yet finished. EOF = ', i_pos
    return
    
end subroutine read_traj
