program main

    implicit none
    integer, parameter:: fid_dat = 101
    character(*), parameter:: file_name = '../../dat/reduced_data.dat'
    integer:: n_traj, i_type, data_n_near_core, ierr_read
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
          data_M_re, data_M_im
    double precision:: w2_max

    n_traj = 0
    w2_max = 0d0
    open( fid_dat, file=trim(file_name) )
    
    do while (1)
        read( fid_dat, '(11(e16.8,1x),i4,1x,i1)', &
              err=100, iostat=ierr_read, end=101 ), &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
              data_M_re, data_M_im, data_n_near_core, i_type;
        n_traj = n_traj + 1
        if( w2_max < (data_M_re*data_M_re + data_M_im*data_M_im) ) then
            w2_max = data_M_re*data_M_re + data_M_im*data_M_im
        end if
        if( mod(n_traj, 1000000) .eq. 0  ) print*, 'line:', n_traj
    end do

100 print*, 'error during reading line', n_traj, ':', ierr_read
101 write(*, '(a,i,1x,a,e16.8)'), 'stop file reading, no of traj = ', n_traj, 'maximum order =', w2_max
    close( fid_dat )
    
end program main
