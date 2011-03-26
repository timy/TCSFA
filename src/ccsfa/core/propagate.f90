!#define PROP_PRNT_INIT_DATA
!#define PROP_PRNT_RK4_DATA

#include '../include/inc_atom.h'
#include '../include/inc_ts_guess.h'


! --------------------------------------------------------------------------------
subroutine propagate_with_single_p0( p0_x, p0_z, ts_guess, &
      ts, amp_M, x0, z0, px_inf, pz_inf, L, n_pass_x, n_pass_z, ierr )


    implicit none;
    double precision, intent(in):: p0_x, p0_z
    double complex, intent(in):: ts_guess
    double complex, external:: solve_ts_from_p0, &
          init_x0, init_z0, init_vx0, init_vz0, &
          action_W_im, action_DDW
    double precision:: t0, tp, x0, z0, vx0, vz0
    double complex:: ts, x0_, z0_, vx0_, vz0_
    double precision, external:: get_tp
    external:: newton_equation_re
    integer:: ierr, n_pass_x, n_pass_z
    double complex:: W_im, W_re, action_W, DDW, amp_M
    double precision:: px_inf, pz_inf, xp, zp, L

    call set_p0( p0_x, p0_z );

    ts = solve_ts_from_p0( ts_guess );    

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

#ifdef PROP_PRNT_INIT_DATA
    print*, 'p0_x', p0_x;
    print*, 'p0_z', p0_z;
    print*, 'ts', ts;
    print*, 'x0', x0;
    print*, 'z0', z0;
    print*, 'vx0', vx0;
    print*, 'vz0', vz0;
#endif


    W_im = action_W_im( ts );
    
    if( dimag(W_im) > 0 ) then
!        ierr = 3;
!        return;
    end if

    DDW = action_DDW( ts );    
    
    call rk4_re( t0, tp, x0, vx0, z0, vz0, &
          newton_equation_re, &
          ierr, W_re, px_inf, pz_inf, &
          xp, zp, L, n_pass_x, n_pass_z);

    action_W = W_im + W_re;

    amp_M = cdexp( -(0d0, 1d0) * action_W ) / DDW;

    ! TEST
!    if(ierr > 0) then
!        print*, 'ierr', ierr
!    end if

#ifdef PROP_PRNT_RK4_DATA
    print*, 'tp', tp
    print*, 'ierr', ierr;
    print*, 'W_re', W_re;
    print*, 'px_inf',px_inf, 'pz_inf',pz_inf;
    print*, 'xp', xp, 'zp', zp
    print*, 'n_pass_x', n_pass_x, 'n_pass_z', n_pass_z
    print*, 'W_im', W_im;
    print*, 'W_re', W_re;
    print*, 'W', action_W;
    print*, 'DDW', DDW;
    print*, 'M', amp_M;
    print*, '----------------------------------------';

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
