program main
    implicit none
    integer:: n_traj, i_type, data_n_near_core, ierr_read, i_z, i_x
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
          data_M_re, data_M_im
    double precision:: dz, dx
    integer, parameter:: nz = 1600
    integer, parameter:: nx = 400
    integer, parameter:: n_type = 4
    double precision:: w(n_type), phi(n_type)
    double precision, parameter:: z0 = -2d0
    double precision, parameter:: z1 = 2d0
    double precision, parameter:: x0 = 0d0
    double precision, parameter:: x1 = 1d0
    double precision:: M_re(nz, nx, n_type), M_im(nz, nx, n_type)
    integer, parameter:: fid_dat = 99
    integer, parameter:: fid_output = 98


    dz = ( z1 - z0 ) / nz
    dx = ( x1 - x0 ) / nx

    M_re = 0d0
    M_im = 0d0
    n_traj = 0

    open( fid_dat, file='../../dat/reduced_data.dat' )
    
    do while(1)
        read( fid_dat, '(11(e16.8,1x),i4,1x,i1)', &
              err=100, iostat=ierr_read, end=101 ), &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
              data_M_re, data_M_im, data_n_near_core, i_type;
        i_z = ceiling( ( data_pz_inf - z0 ) / dz )
        i_x = ceiling( ( data_px_inf - x0 ) / dx )
        if( i_z .lt. 1 .or. i_z .gt. nz ) cycle
        if( i_x .lt. 1 .or. i_x .gt. nx ) cycle
        !print*, 'i_z',i_z,'i_x',i_x, data_pz_inf, data_px_inf
        M_re(i_z,i_x,i_type) = M_re(i_z,i_x,i_type) + data_M_re
        M_im(i_z,i_x,i_type) = M_im(i_z,i_x,i_type) + data_M_im
        n_traj = n_traj + 1
        if( mod(n_traj,100000) .eq. 0 ) then
            print*, "count = ", n_traj
        end if
                
    end do

100 print*, 'error during reading line', n_traj, ':', ierr_read
101 write(*, '(a,i,1x,a,e16.8)'), 'stop file reading, no of traj = ', n_traj

    close( fid_dat )
    
    open( fid_output, file='../../dat/phase.dat' )
    do i_x = 1, nx
        do i_z = 1, nz

            do i_type = 1, n_type
                w(i_type) = M_re(i_z,i_x,i_type)*M_re(i_z,i_x,i_type) &
                      + M_im(i_z,i_x,i_type)*M_im(i_z,i_x,i_type)
                phi(i_type) = datan2( M_im(i_z,i_x,i_type), M_re(i_z,i_x,i_type) )
            end do
            write( fid_output, '(9(e16.8,1x))' ), &
                  sum( M_re(i_z,i_x,:) )**2 + sum( M_im(i_z,i_x,:) )**2, &
                  (w(i_type), i_type=1,n_type), &
                  (phi(i_type),i_type=1,n_type)
        end do
    end do

    close( fid_output )

    write(*,'(a)'), "hasta la vista"
end program main
