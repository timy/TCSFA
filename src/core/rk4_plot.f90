#include '../include/inc_misc.h'
#ifdef MISC_PLOT
#include '../include/inc_plot_rk4.h'

module mod_rk4_plot

    implicit none

    integer:: os
    character(len=64):: file_name

end module mod_rk4_plot

! ----------------------------------------------------------------------
subroutine rk4_plot_init(i)
    use mod_rk4_plot, only: os, file_name
    implicit none
    integer, intent(in):: i

    os = RK4_TRAJ_PLOT_FILE_ID + i - 1
    write( file_name, '(a,i5,a)' ), RK4_TRAJ_PLOT_FILE_NAME, i, '.dat'
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
subroutine rk4_plot_write( n_step, t, y, h, n_substep )
    use mod_rk4_plot, only: os
    implicit none
    double precision:: t, y(4), h
    integer, intent(in):: n_step, n_substep
    ! order of z, x, vz, vx
    write( os, '(i6,x,6(e16.8,x),i6)' ), n_step, t, y(3), y(1), y(4), y(2), h, n_substep
    return
end subroutine rk4_plot_write

! ----------------------------------------------------------------------
subroutine rk4_plot_close_file()
    use mod_rk4_plot, only: os
    implicit none
    close( os )
    return
end subroutine rk4_plot_close_file

#endif
