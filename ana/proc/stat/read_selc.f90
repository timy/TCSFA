subroutine read_selc( n_traj, fid_selc )
    implicit none
    integer, intent(in):: n_traj, fid_selc
    integer, parameter:: n_bin = 100
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
          data_M_re, data_M_im
    integer:: data_n_near_core, i_type, i_pos, os
    double precision, allocatable:: ts_re(:), ts_im(:), Mp_re(:), Mp_im(:), &
          x0(:), z0(:), p0_x(:), p0_z(:), L(:), w(:)
    integer, allocatable:: n_near_core(:), types(:)
    character(*), parameter:: dir_dat = "../../dat/"

    allocate( ts_re(n_traj) ); allocate( ts_im(n_traj) )
    allocate( Mp_re(n_traj) ); allocate( Mp_im(n_traj) )
    allocate( x0(n_traj) ); allocate( z0(n_traj) )
    allocate( p0_x(n_traj) ); allocate( p0_z(n_traj) )
    allocate( n_near_core(n_traj) )
    allocate( L(n_traj) ); allocate( types(n_traj) )
    allocate( w(n_traj) )

    do i_pos = 1, n_traj
        read( fid_selc, '(11(e16.8e3,1x),i4,1x,i1)' ), &
              data_px_0, data_pz_0, data_ts_re, data_ts_im, &
              data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
              data_M_re, data_M_im, data_n_near_core, i_type

        ts_re(i_pos) = data_ts_re; ts_im(i_pos) = data_ts_im
        Mp_re(i_pos) = data_M_re; Mp_im(i_pos) = data_M_im
        x0(i_pos) = data_x_0; z0(i_pos) = data_z_0
        p0_x(i_pos) = data_px_0; p0_z(i_pos) = data_pz_0
        n_near_core(i_pos) = data_n_near_core
        L(i_pos) = data_L; types(i_pos) = i_type
        w(i_pos) = dsqrt( Mp_re(i_pos)**2 + Mp_im(i_pos)**2 )
    end do

    os = 110
    open( os, file=dir_dat//'hist_ts_re.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_ts_re.dat ...'
    call histogram_w( n_traj, ts_re, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_ts_im.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_ts_im.dat ...'
    call histogram_w( n_traj, ts_im, w, types, n_bin, os )
    close( os )
    
    os = 110
    open( os, file=dir_dat//'hist_Mp_re.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_Mp_re.dat ...'
    call histogram_w( n_traj, Mp_re, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_Mp_im.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_Mp_im.dat ...'
    call histogram_w( n_traj, Mp_im, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_x0.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_x0.dat ...'
    call histogram_w( n_traj, x0, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_z0.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_z0.dat ...'
    call histogram_w( n_traj, z0, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_p0x.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_p0x.dat ...'
    call histogram_w( n_traj, p0_x, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_p0z.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_p0z.dat ...'
    call histogram_w( n_traj, p0_z, w, types, n_bin, os )
    close( os )

    os = 110
    open( os, file=dir_dat//'hist_L.dat' )
    write(*, '(a)'), 'writing to '//dir_dat//'hist_L.dat ...'
    call histogram_w( n_traj, L, w, types, n_bin, os )
    close( os )

    deallocate( ts_re ); deallocate( ts_im )
    deallocate( Mp_re ); deallocate( Mp_im )
    deallocate( x0 ); deallocate( z0 ); 
    deallocate( p0_x ); deallocate( p0_z ); 
    deallocate( n_near_core )
    deallocate( L ); deallocate( types )
    deallocate( w )
    
    return
end subroutine read_selc

subroutine histogram( n, a, n_bin, os )
    implicit none
    integer, intent(in):: n, n_bin
    double precision, intent(in):: a(n)
    double precision:: bin(n_bin), a_min, a_max, a_lower, da
    integer:: i, i_a, os
    
    a_min = minval(a)
    a_max = maxval(a)
    write(*,'(2(a, e16.8e3, 2x))'), 'min = ', a_min, 'max = ', a_max
    da = ( a_max - a_min ) / ( n_bin - 1d0 )
    if( da .eq. 0d0 ) then
        write(*, '(a)'), 'histogram: no interval, neglected'
        return
    end if
    a_lower = a_min - 0.5 * da
    bin = 0

    do i = 1, n
        i_a = ceiling( ( a(i) - a_lower ) / da )
        if( i_a < 0 .or. i_a > n_bin ) then
            print*, 'i_a = ', i_a
            print*, "a_min = ", a_min, "a_max = ", a_max
            print*, "a(i):", i, a(i)
        end if
        bin(i_a) = bin(i_a) + 1.0
    end do

    print*, "here"

    if( os > 0 ) then
        do i = 1, n_bin
            write( os, '(2(e16.8e3,1x))' ), a_min + (i-1)*da, bin(i)
        end do
    end if

    return
end subroutine histogram

subroutine histogram_w( n, a, w, types, n_bin, os )
    implicit none
    integer:: n, n_bin
    double precision, intent(in):: a(n), w(n)
    integer:: types(n)
    integer, parameter:: n_type = 5
    double precision:: bin(n_bin, n_type), a_min, a_max, a_lower, da
    integer:: i, i_a, os, j
    
    a_min = minval(a)
    a_max = maxval(a)
    write(*,'(2(a, e16.8e3, 2x))'), 'min = ', a_min, 'max = ', a_max
    da = ( a_max - a_min ) / ( n_bin - 1d0 )
    if( da .eq. 0d0 ) then
        write(*, '(a)'), 'histogram: no interval, neglected'
        return
    end if
    a_lower = a_min - 0.5d0 * da
    bin = 0d0

    do i = 1, n
        i_a = ceiling( ( a(i) - a_lower ) / da )
        bin(i_a, types(i)) = bin(i_a, types(i)) + w(i)
    end do
        
    if( os > 0 ) then
        do i = 1, n_bin
            write( os, '(5(e16.8e3,1x))' ), a_min + (i-1)*da, ( bin(i,j), j=1,4 )
        end do
    end if

    return
end subroutine histogram_w

