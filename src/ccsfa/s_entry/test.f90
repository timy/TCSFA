#include "../include/inc_field.h"

program main

    use mod_pulse, only: set_pulse, set_pulse_t0
    
    implicit none

    double precision:: p0_x, p0_z, pinf_x, pinf_z
    double complex:: ts_guess, W_im, W_re
    double complex, external:: pulse_A_z, pulse_A_x, pulse_E_z, pulse_E_x, alpha_x, alpha_z

    call set_pulse( E0, OM, NC, XI, PH );
    call set_pulse_t0( (238.769d0, 0d0) );
    call plot_pulse

!    call read_data_and_plot_traj()

    p0_x = -0.42901813d0;
    p0_z = -3.0137534d0;
    ts_guess = ( 470.60866d0, 120.860932d0 );

    call disp_one_traj( p0_x, p0_z, ts_guess )
    
end program main







subroutine disp_one_traj( p0_x, p0_z, ts_guess )

    implicit none
    double precision, intent(in):: p0_x, p0_z
    double complex, intent(in):: ts_guess
    double complex:: ts, amp_M
    double precision:: x0, z0, px_inf, pz_inf
    integer:: n_pass_x, n_pass_z, ierr
    integer:: fid = 6


    write(fid, '(a)'), 'Initial Conditions'
    print*, 'p0_x', p0_x, 'p0_z', p0_z
    print*, 'ts_guess:', ts_guess
    
    call propagate_with_single_p0( &
          p0_x, p0_z, ts_guess, &
          ts, amp_M, x0, z0, px_inf, pz_inf, &
          n_pass_x, n_pass_z, ierr )    

!!$    write(fid, '(a)'), 'Results'
!!$    print*, 'ts', ts;
!!$    print*, 'amp_M', amp_M
!!$    print*, 'x0', x0, 'z0', z0
!!$    print*, 'px', px_inf, 'pz', pz_inf
!!$    print*, 'n_pass_x', n_pass_x, 'n_pass_z', n_pass_z
!!$    print*, 'ierr', ierr;

    return;
    
end subroutine disp_one_traj
