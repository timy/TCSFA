#include '../../../src/include/inc_field.h'
#include '../../../src/include/inc_grid_boundary.h'
#include '../../../src/include/inc_grid_size.h'

! ---------------------------------------------------------------------- !
! "main.trajs.f90" is used to generate trajectory(ies) from "selc_*.dat" !
! ---------------------------------------------------------------------- !

program main
    implicit none
    integer, parameter:: file_index = 1
    integer, parameter:: n_traj = 100
    integer, parameter:: fid_selc = 20
    integer, parameter:: fid_info = 21
    double precision:: x0, z0, px_inf, pz_inf, L, w, w_max
    double complex:: ts, Mp, Mp_total
    integer:: n_pass_x, n_pass_z, n_step, ierr
    ! data read from selc_*.dat
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
          data_M_re, data_M_im
    integer:: count, i_pos, i_type, ierr_read, data_n_step
    character(len=128):: filename

    ! read initial data from file "selc_*.dat"
    write( filename, '(a,i3,a)' ), '../../dat/selc/selc_', file_index, '.dat'
    open( fid_selc, file=filename )
    open( fid_info, file='../../dat/traj/info.dat' )
    count = 0
    i_pos = 0
    do while( count < n_traj )
       ! read in raw data
       read( fid_selc, '(11(e16.8,1x),i4,1x,i1)', &
            err=101, iostat=ierr_read, end=105 ), &
            data_px_0, data_pz_0, data_ts_re, data_ts_im, &
            data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
            data_M_re, data_M_im, data_n_step, i_type
       
       i_pos = i_pos + 1

       w = data_M_re**2 + data_M_im**2
       if( w_max < w ) w_max = w

       ! --------------------------------------------------------------------
       ! set the filter, to only plot trajectores meeting the criteria
       !if( i_type .ne. 2 ) cycle
       if( w .le. 2d-5 ) cycle
       ! --------------------------------------------------------------------
       count = count + 1
       Mp_total = Mp_total + Mp

       call set_p0( data_px_0, data_pz_0 )  
       call propagate_with_single_p0( data_px_0, data_pz_0, &
            dcmplx(data_ts_re, data_ts_im), &
            ts, Mp, x0, z0, px_inf, pz_inf, L, &
            n_pass_x, n_pass_z, n_step, ierr, count )
       write(*, '(a,i3)'), 'traj: ', count

       write( fid_info, '(i4,1x,e15.8,1x,i1)' ), count, &
             data_M_re**2+data_M_im**2, i_type

       cycle
101    print*, 'error during reading line', i_pos, ':', ierr_read;
    end do
    ! close files
    go to 110

105 write(*, *), 'EOF = ', i_pos
    
110 close( fid_selc )
    close( fid_info )

    write(*, '(a, e15.8, 1x, e15.8)'), "total Mp = ", Mp_total
    write(*, '(a, e15.8)'),  "w_max = ", w_max
    write(*, '(a, i4)'), 'output traj no. :', count
    write(*, '(a)'), "hasta la vista"
end program main
