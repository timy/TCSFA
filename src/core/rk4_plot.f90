#include '../include/inc_misc.h'
#ifdef MISC_PLOT_TRAJ
#include '../include/inc_plot_rk4.h'

module mod_rk4_plot
    implicit none
    integer:: os
    integer:: flag
    character(len=256):: file_name
end module mod_rk4_plot

! ----------------------------------------------------------------------
subroutine rk4_plot_init(flag_, i)
    use mod_rk4_plot, only: os, flag, file_name
    implicit none
    integer, intent(in):: flag_, i
    ! flag indicate the plot type, i assign a ID for the traj
    flag = flag_
    os = RK4_TRAJ_PLOT_FILE_ID + i - 1
    if( flag .eq. 0 ) then ! real trajectory
        write( file_name, '(a,a,i5,a)' ), RK4_TRAJ_PLOT_FILE_NAME, '_re_', i, '.dat'
        write( *, '(a, a)' ),  're_traj is written to ', file_name
    elseif( flag .eq. 1 ) then ! complex trajectory
        write( file_name, '(a,a,i5,a)' ), RK4_TRAJ_PLOT_FILE_NAME, '_im_', i, '.dat'
        write( *, '(a, a)' ),  'im_traj is written to ', file_name        
    end if
    return
end subroutine rk4_plot_init

! ----------------------------------------------------------------------
subroutine rk4_plot_open_file()
    use mod_rk4_plot, only: os, file_name
    implicit none
    open( os, file = trim( file_name ) )    
    return
end subroutine rk4_plot_open_file

! ----------------------------------------------------------------------
subroutine rk4_plot_write_real( n_step, t, y, h, n_substep )
    use mod_rk4_plot, only: os, flag
    implicit none
    double precision, intent(in):: t, h
    double precision, intent(in):: y(4)
    integer, intent(in):: n_step, n_substep
    ! order of z, x, vz, vx
    write( os, '(i6,x,6(e16.8,x),i6)' ), n_step, t, y(3), y(1), y(4), y(2), h, n_substep
    return
end subroutine rk4_plot_write_real

! ----------------------------------------------------------------------
subroutine rk4_plot_write_cmplx( n_step, t, y, h, n_substep )
    use mod_rk4_plot, only: os, flag
    implicit none
    double precision, intent(in):: t, h
    double complex, intent(in):: y(4)
    integer, intent(in):: n_step, n_substep
    ! order of z, x, vz, vx
    write( os, '(i6,x,10(e16.8,x),i6)' ), n_step, t, y(3), y(1), y(4), y(2), h, n_substep
    return
end subroutine rk4_plot_write_cmplx

! ----------------------------------------------------------------------
subroutine rk4_plot_close_file()
    use mod_rk4_plot, only: os
    implicit none
    close( os )
    return
end subroutine rk4_plot_close_file

#endif
