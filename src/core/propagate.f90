!#define PROP_PRNT_INIT_DATA
!#define PROP_PRNT_RK4_DATA

#include '../include/inc_atom.h'
#include '../include/inc_ts_guess.h'
#include '../include/inc_misc.h'
#include '../include/inc_field.h'


#define SUB_PTB_PROP sub_ptb_prop
!#define SUB_PTB_PROP sub_ptb_prop_0

! --------------------------------------------------------------------------------
subroutine propagate_with_single_p0( p0_x, p0_z, ts_guess, &
      ts, amp_M, x0, z0, px_inf, pz_inf, L, n_near_core, err_spe, ierr, tag )


    implicit none;
    double precision, intent(in):: p0_x, p0_z
    double complex, intent(in):: ts_guess
    integer, intent(in):: tag
    double complex, external:: solve_ts_from_p0, &
          init_x0, init_z0, init_vx0, init_vz0, &
          action_DDW!, action_W_sub_0_traj_0, action_W_sub_r_rcpr, action_W_sub_ptb
    double precision:: t0, x0, z0, vx0, vz0
    double complex:: ts, x0_, z0_, vx0_, vz0_
    external:: newton_equation_re
    integer:: ierr, n_near_core
    double complex:: W_sub, W_re, action_W, DDW, amp_M !W_sub_r_rcpr, 
    double precision:: px_inf, pz_inf, L
    double precision:: err_spe
!#if MISC_PRINT > 2
    double complex, external:: PULSE_E_Z, PULSE_E_X, PULSE_A_Z, PULSE_A_X
!#endif

    call set_p0( p0_x, p0_z );
    ts = solve_ts_from_p0( ts_guess, ierr );    
    if( ierr > 0 ) then
        return
    end if

    call SUB_PTB_PROP( ts, ierr, z0_, x0_, vz0_, vx0_, w_sub, err_spe, tag )
    if( ierr > 0 ) return


! --------------------------------------------------------------------------------
! W_sub and sub-barrier trajectory
#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to calculate action for sub-barrier '
#endif

!    W_sub = action_W_sub_0_traj_0(ts)

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_sub                 (', W_sub,                    ' )'
#endif

    ! add the contribution from 1/r 
!    W_sub_r_rcpr = action_W_sub_r_rcpr(ts)

!#if MISC_PRINT > 2
!    write(*,'(a)'), ''
!    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_sub_r_rcpr      (', W_sub_r_rcpr,         '  )'
!#endif

!   W_sub = W_sub + W_sub_r_rcpr
!    W_sub = action_W_sub_ptb(ts)
! --------------------------------------------------------------------------------
! DDW
#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to calculate DDW"(ts)'
#endif

    DDW = action_DDW( ts );    

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'DDW       (', DDW,         '  )'
#endif

! --------------------------------------------------------------------------------
! initial conditions for real propagation
#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to obtain boundary values from ts '
#endif

    ! x0_ = init_x0( ts );
    ! z0_ = init_z0( ts );
    ! vx0_ = init_vx0( ts );
    ! vz0_ = init_vz0( ts );

    t0 = dreal( ts );

    x0 = dreal( x0_ );
    z0 = dreal( z0_ );
    vx0 = dreal( vx0_ );
    vz0 = dreal( vz0_ );

#if MISC_PRINT > 2

    write(*,'(a)'), ''
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'ts        (', ts,           '  )'
    write(*,'(4x, a, f15.8)'),                'tp        ', Tp
    write(*,'(2(4x, a, f15.8))'),             'x(t0)     ', x0,       'z(t0)     ', z0
    write(*,'(2(4x, a, f15.8))'),             'vx(t0)    ', vx0,       'vz(t0)    ', vz0
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ez(ts)    (', PULSE_E_Z( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ex(ts)    (', PULSE_E_X( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Az(ts)    (', PULSE_A_Z( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ax(ts)    (', PULSE_A_X( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ez(t0)    (', PULSE_E_Z( dcmplx(t0,0d0) ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ex(t0)    (', PULSE_E_X( dcmplx(t0,0d0) ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Az(t0)    (', PULSE_A_Z( dcmplx(t0,0d0) ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ax(t0)    (', PULSE_A_X( dcmplx(t0,0d0) ),          '  )'

#endif

! --------------------------------------------------------------------------------
#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to calculate real trajectory with rk4'
#endif

    call rk4_prop( t0, Tp, x0, vx0, z0, vz0, ierr, W_re, px_inf, pz_inf, L, n_near_core, tag );

    if( ierr == 0 ) then
        action_W = W_sub + W_re;
        ! this one is for 1s state
        amp_M = cdexp( dcmplx(0d0, 1d0) * action_W ) / DDW;
        ! this one is for 2p state
!        amp_M = ((PULSE_A_z(ts)+p0_z)/dsqrt(2d0*IONIZATION_IP)) * cdexp( -dcmplx(0d0, 1d0) * action_W ) / DDW
 
! if you want to check the sub-barrier trajectory of the 1st order
!        call plot_sub_traj_ptb_1(ts)
! end sub-barrier trajectory demo
    else
        amp_M = dcmplx(0d0, 0d0);
    end if
        

#if MISC_PRINT > 2

    write(*,'(a)'), ''
    write(*,'(4x, a, 9x, i6)'),               'ierr      ', ierr;
    write(*,'(2(4x, a, f15.8))'),             'px_inf    ', px_inf,   'pz_inf    ', pz_inf
!    write(*,'(2(4x, a, f15.8))'),             'x(tp)     ', xp,       'z(tp)     ', zp
!    write(*,'(2(4x, a, 9x, i6))'),            'n_pass_x  ', n_pass_x, 'n_pass_z  ', n_pass_z
    write(*,'(4x, a, f15.8)'),                'L         ', L
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_sub     (', W_sub,        '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_re       (', W_re,          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W            (', action_W,    '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'DDW       (', DDW,          '  )'
    write(*,'(4x, a, e15.8, 2x, e15.8, a)'),  'Mp         (', amp_M,        '  )'

#endif

    return;
end subroutine propagate_with_single_p0



subroutine newton_equation_re( ne, t, y, dy )
    
    implicit none;
    integer, intent(in):: ne;
    double precision, intent(in):: t, y(ne);
    double precision, intent(out):: dy(ne);
    double complex, external:: PULSE_E_Z, PULSE_E_X;
    double precision:: r, r3, charge;

    !double precision, parameter:: alpha = 5.0
    !double precision:: effective_pot_term;
    ! y(1) = x; y(2) = vx; y(3) = z; y(4) = vz;

    charge = ATOM_CHARGE_Z;
    
    r = dsqrt( y(1)*y(1) + y(3)*y(3) );
    r3 = r * r * r;
    !effective_pot_term = -( 1d0 + dexp(-r*alpha) * ( charge - 1d0 ) * ( 1d0 + alpha*r ) ) / r3
    dy(1) = y(2); 
    dy(2) = - dreal( PULSE_E_X( dcmplx(t) ) ) - y(1) * charge / r3;
    !dy(2) = -dreal( PULSE_E_X( dcmplx(t) ) ) + y(1) * effective_pot_term
    dy(3) = y(4);
    dy(4) = - dreal( PULSE_E_Z( dcmplx(t) ) ) - y(3) * charge / r3;
    !dy(4) = -dreal( PULSE_E_Z( dcmplx(t) ) ) + y(3) * effective_pot_term

    return;
end subroutine Newton_equation_re

! ! ////////////////////////////////////////////////////////////////////////////////
! ! x(t_0)
! double complex function init_x0( ts )

!     implicit none
!     double complex, intent(in):: ts
!     double complex, external:: sub_traj_x_0, sub_traj_x_1

!     init_x0 = sub_traj_x_0( dcmplx( dreal(ts), 0d0 ), ts ) + sub_traj_x_1( 0d0 )
!     return;
! end function init_x0

! ! ////////////////////////////////////////////////////////////////////////////////
! ! z(t_0)
! double complex function init_z0( ts )

!     implicit none;
!     double complex, intent(in):: ts;
!     double complex, external:: sub_traj_z_0!, sub_traj_z_1

!     init_z0 = sub_traj_z_0( dcmplx( dreal(ts), 0d0 ), ts )! + sub_traj_z_1( 0d0 )
!     return;
! end function init_z0

! ! ////////////////////////////////////////////////////////////////////////////////
! ! vx(t_0)
! double complex function init_vx0( ts )

!     implicit none;
!     double complex, intent(in):: ts;
!     double complex, external:: sub_traj_vx_0, sub_traj_vx_1;

!     init_vx0 = sub_traj_vx_0( dcmplx(dreal(ts), 0d0) ) !+ sub_traj_vx_1( 0d0, ts )
!     return;
! end function init_vx0

! ! ////////////////////////////////////////////////////////////////////////////////
! ! vz(t_0)
! double complex function init_vz0( ts )
!     use mod_p0;

!     implicit none;
!     double complex, intent(in):: ts;
!     double complex, external:: sub_traj_vz_0

!     init_vz0 = sub_traj_vz_0( dcmplx(dreal(ts), 0d0) ) !+ sub_traj_vz_1( 0d0, ts )
!     return;
! end function init_vz0
