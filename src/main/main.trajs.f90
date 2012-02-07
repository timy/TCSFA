#include '../include/inc_field.h'
#include '../include/inc_ts_guess.h'
#include '../include/inc_grid_boundary.h'
#include '../include/inc_grid_size.h'

! ---------------------------------------------------------------------- !
! "main.trajs.f90" is used to generate trajectory(ies) from "selc_*.dat" !
! ---------------------------------------------------------------------- !

program main
    use mod_pulse
    implicit none
    integer, parameter:: file_index = 40
    integer, parameter:: n_traj = 2
    character(len=40):: filename
    ! data read from "selc_*.dat"
    double precision:: data_px_0, data_pz_0, data_ts_re, data_ts_im, &
          data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
          data_M_re, data_M_im
    integer:: data_n_near_core, i_type
    ! data obtained from re-calculation
    double precision:: x0, z0, px_inf, pz_inf, L
    double complex:: ts, Mp
    integer:: n_pass_x, n_pass_z, n_near_core, ierr
    integer, parameter:: fid_selc = 20
    integer:: i_pos, ierr_read, count

    ! print the laser parameters
    print*, 'E0', E0;
    print*, 'OM', OM;
    print*, 'NC', NC;
    print*, 'XI', XI;
    print*, 'PH', PH;

    call set_pulse( E0, OM, NC, XI, PH )
    call set_pulse_t0( ( 1d0, 0d0 ) )
    call plot_pulse()

    ! read initial data from file "selc_*.dat"
    write( filename, '(a,i3,a)' ), 'ana/dat/selc_', file_index, '.dat'
    open( fid_selc, file=filename )
    
    count = 0
    i_pos = 0
    do while( count < n_traj )
       ! read in raw data
       read( fid_selc, '(11(e16.8,1x),i4,1x,i1)', &
            err=101, iostat=ierr_read, end=105 ), &
            data_px_0, data_pz_0, data_ts_re, data_ts_im, &
            data_x_0, data_z_0, data_px_inf, data_pz_inf, data_L, &
            data_M_re, data_M_im, data_n_near_core, i_type
       
       i_pos = i_pos + 1
       ! set the filter, to only plot trajectores meeting the criteria
       if( i_type .ne. 3 ) cycle
       
       count = count + 1

       call set_p0( data_px_0, data_pz_0 )  
       call propagate_with_single_p0( data_px_0, data_pz_0, &
            dcmplx(data_ts_re, data_ts_im), &
            ts, Mp, x0, z0, px_inf, pz_inf, L, &
            n_pass_x, n_pass_z, n_near_core, ierr, count )
       write(*, '(a,i3)'), 'traj: ', count
       cycle
101    print*, 'error during reading line', i_pos, ':', ierr_read;
    end do
    ! close files
    close( fid_selc )
    stop "hasta la vista"

105 write(*, *), 'EOF = ', i_pos
    stop ":("
  ! stop here!!!

end program main
