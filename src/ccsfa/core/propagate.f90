!#define PROP_PRNT_INIT_DATA
!#define PROP_PRNT_RK4_DATA

#include '../include/inc_atom.h'
#include '../include/inc_ts_guess.h'
#include '../include/inc_misc.h'

! --------------------------------------------------------------------------------
subroutine propagate_with_single_p0( p0_x, p0_z, ts_guess, &
      ts, amp_M, x0, z0, px_inf, pz_inf, L, n_pass_x, n_pass_z, ierr )


    implicit none;
    double precision, intent(in):: p0_x, p0_z
    double complex, intent(in):: ts_guess
    double complex, external:: solve_ts_from_p0, &
          init_x0, init_z0, init_vx0, init_vz0, &
          action_W_im, action_W_im_num, action_DDW
    double precision:: t0, tp, x0, z0, vx0, vz0
    double complex:: ts, x0_, z0_, vx0_, vz0_
    double precision, external:: get_tp
    external:: newton_equation_re
    integer:: ierr, n_pass_x, n_pass_z
    double complex:: W_im, W_re, action_W, DDW, amp_M
    double precision:: px_inf, pz_inf, xp, zp, L
#if MISC_PRINT > 2
    double complex, external:: pulse_E_z, pulse_E_x, pulse_A_z, pulse_A_x
#endif


#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(a)'), repeat('#', 50)
    write(*,'(a)'), 'propagate.f90: propagate_with_single_p0 '
#endif

    call set_p0( p0_x, p0_z );

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to solve ts '
#endif

    ts = solve_ts_from_p0( ts_guess );    

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to obtain boundary values from ts '
#endif

    x0_ = init_x0( ts );
    z0_ = init_z0( ts );
    vx0_ = init_vx0( ts );
    vz0_ = init_vz0( ts );

    t0 = dreal( ts );
    tp = get_tp();

    x0 = dreal( x0_ );
    z0 = dreal( z0_ );
    vx0 = dreal( vx0_ );
    vz0 = dreal( vz0_ );

#if MISC_PRINT > 2

    write(*,'(a)'), ''
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'ts        (', ts,           '  )'
    write(*,'(4x, a, f15.8)'),                'tp        ', tp
    write(*,'(2(4x, a, f15.8))'),             'x(t0)     ', x0,       'z(t0)     ', z0
    write(*,'(2(4x, a, f15.8))'),             'vx(t0)    ', vx0,       'vz(t0)    ', vz0
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ez(ts)    (', pulse_E_z( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ex(ts)    (', pulse_E_x( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Az(ts)    (', pulse_A_z( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ax(ts)    (', pulse_A_x( ts ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ez(t0)    (', pulse_E_z( dcmplx(t0,0d0) ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ex(t0)    (', pulse_E_x( dcmplx(t0,0d0) ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Az(t0)    (', pulse_A_z( dcmplx(t0,0d0) ),          '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'Ax(t0)    (', pulse_A_x( dcmplx(t0,0d0) ),          '  )'

#endif

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to calculate action for sub-barrier '
#endif

    W_im = action_W_im( ts );
    
#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_im      (', W_im,         '  )'
#endif

    if( dimag(W_im) > 0 ) then
!        ierr = 3;
!        return;
    end if

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to calculate W"(ts)'
#endif

    DDW = action_DDW( ts );    
    
#if MISC_PRINT > 2

    write(*,'(a)'), ''
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'DDW       (', DDW,         '  )'

#endif

#if MISC_PRINT > 2
    write(*,'(a)'), ''
    write(*,'(2x, a)'), '* start to calculate real trajectory with rk4'
#endif

    call rk4_re( t0, tp, x0, vx0, z0, vz0, &
          newton_equation_re, &
          ierr, W_re, px_inf, pz_inf, &
          xp, zp, L, n_pass_x, n_pass_z);

    action_W = W_im + W_re;

    amp_M = cdexp( -(0d0, 1d0) * action_W ) / DDW;


#if MISC_PRINT > 2

    write(*,'(a)'), ''
    write(*,'(4x, a, 9x, i6)'),               'ierr      ', ierr;
    write(*,'(2(4x, a, f15.8))'),             'px_inf    ', px_inf,   'pz_inf    ', pz_inf
    write(*,'(2(4x, a, f15.8))'),             'x(tp)     ', xp,       'z(tp)     ', zp
    write(*,'(2(4x, a, 9x, i6))'),            'n_pass_x  ', n_pass_x, 'n_pass_z  ', n_pass_z
    write(*,'(4x, a, f15.8)'),                'L         ', L
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_im     (', W_im,        '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W_re     (', W_re,        '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'W        (', action_W,    '  )'
    write(*,'(4x, a, f15.8, 2x, f15.8, a)'),  'DDW      (', DDW,         '  )'
    write(*,'(4x, a, e15.8, 2x, e15.8, a)'),  'Mp       (', amp_M,       '  )'

#endif

    return;
end subroutine propagate_with_single_p0




















subroutine newton_equation_re( ne, t, y, dy )
    
    implicit none;
    integer, intent(in):: ne;
    double precision, intent(in):: t, y(ne);
    double precision, intent(out):: dy(ne);
    double complex, external:: pulse_E_z, pulse_E_x;
    double precision:: r, r3, charge;
    ! y(1) = x; y(2) = vx; y(3) = z; y(4) = vz;


    charge = ATOM_CHARGE_Z;
    r = dsqrt( y(1)*y(1) + y(3)*y(3) );
    r3 = r * r * r;
    dy(1) = y(2);
    dy(2) = - dreal( pulse_E_x( dcmplx(t) ) ) - y(1) * charge / r3;
    dy(3) = y(4);
    dy(4) = - dreal( pulse_E_z( dcmplx(t) ) ) - y(3) * charge / r3;
    return;
end subroutine Newton_equation_re








double complex function init_x0( ts )

    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: alpha_x;
    double complex:: t0;

    t0 = dcmplx( dreal( ts ), 0d0 );
    init_x0 = alpha_x( t0 ) - dreal( alpha_x( ts ) );

    return;
end function init_x0


double complex function init_z0( ts )

    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: alpha_z;
    double complex:: t0;

    t0 = dcmplx( dreal( ts ), 0d0 );
    init_z0 = alpha_z( t0 ) - dreal( alpha_z( ts ) );

    return;
end function init_z0


double complex function init_vx0( ts )
    use mod_p0;

    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: pulse_A_x;
    double complex:: t0;

    t0 = dcmplx( dreal( ts ), 0d0 );
    init_vx0 = p0_x + pulse_A_x( t0 );

    return;
end function init_vx0


double complex function init_vz0( ts )
    use mod_p0;

    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: pulse_A_z;
    double complex:: t0;

    t0 = dcmplx( dreal( ts ), 0d0 );
    init_vz0 = p0_z + pulse_A_z( t0 );

    return;
end function init_vz0
