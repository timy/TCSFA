#include '../../../include/inc_atom.h'
#include '../../../include/inc_field.h'

! the sub-barrier action for the 0th-order trajectory without -1/r term
! for the linearly polarized field along the z-axis
subroutine action_W_sub_0_0_0_ana( ts, w_0, w_r_recp )
    use mod_p0, only: p0_x, p0_z
    implicit none
    double complex, intent(in):: ts
    double complex, external:: PULSE_ALPHA_Z, PULSE_ALPHA2_Z
    double complex:: t0, w_0, w_r_recp

    t0 = dcmplx( dreal(ts), 0d0 )
    w_0 = &
         (0.5d0 * (p0_z*p0_z+p0_x*p0_x) + IONIZATION_IP ) * ( ts - t0 ) &
         + p0_z * ( PULSE_ALPHA_Z(ts) - PULSE_ALPHA_Z(t0) ) &
         + 0.5d0 * ( PULSE_ALPHA2_Z(ts) - PULSE_ALPHA2_Z(t0) )
    w_r_recp = dcmplx( 0d0, 0d0 )
    return;
end subroutine  action_W_sub_0_0_0_ana

! the sub-barrier action for the 0th-order trajectory with -1/r term
! for the linearly polarized field along the z-axis
subroutine action_w_sub_1_0_0_ana(ts, w_0, w_r_recp, w_r_recp_abs)
    use mod_p0, only: p0_x, p0_z
    implicit none
    double complex, intent(in):: ts
    double complex, external:: PULSE_ALPHA_Z, PULSE_ALPHA2_Z
    double complex, external:: sub_int_r_reciprocal, sub_int_r_reciprocal_abs
    double complex, external:: simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: t0, w_0, w_r_recp, w_r_recp_abs
    double precision:: ti

    t0 = dcmplx( dreal(ts), 0d0 )
    ti = dimag(ts)
    w_0 = &
         (0.5d0 * (p0_z*p0_z+p0_x*p0_x) + IONIZATION_IP ) * ( ts - t0 ) &
         + p0_z * ( PULSE_ALPHA_Z(ts) - PULSE_ALPHA_Z(t0) ) &
         + 0.5d0 * ( PULSE_ALPHA2_Z(ts) - PULSE_ALPHA2_Z(t0) )

    w_r_recp = - eye * simpson_sub(ti, 0d0, 200, ts, sub_int_r_reciprocal)
    w_r_recp_abs = - eye * simpson_sub(ti, 0d0, 200, ts, sub_int_r_reciprocal_abs)
    return
end subroutine  action_w_sub_1_0_0_ana

! Integrand - Z/r
double complex function sub_int_r_reciprocal( t, ts ) result(w_r_recp)
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex:: tc, r, z0, x0
    double complex, external:: sub_traj_z_0, sub_traj_x_0
    
    tc = dcmplx( dreal(ts), t )
    z0 = sub_traj_z_0( tc, ts )
    x0 = sub_traj_x_0( tc, ts )
    
    r = (z0*z0 + x0*x0)**0.5
    w_r_recp = - ATOM_CHARGE_Z / r
!    w_r_recp = - 1d0 / r
        
    return
end function sub_int_r_reciprocal

! Integrand - Z/abs(r)
double complex function sub_int_r_reciprocal_abs( t, ts ) result(w_r_recp_abs)
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex:: r, tc, r_abs, z0, x0
    double complex, external:: sub_traj_z_0, sub_traj_x_0
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double precision, parameter:: PI = 2d0 * dasin(1d0)
    integer, parameter:: k = 2
    double precision:: phi, x, y

    tc = dcmplx( dreal(ts), t )
    z0 = sub_traj_z_0( tc, ts )
    x0 = sub_traj_x_0( tc, ts )
    
    r = (z0*z0 + x0*x0)**0.5
!    r_abs = cdabs(r)
    w_r_recp_abs = - dreal( ATOM_CHARGE_Z / r )
!    w_r_recp_abs = - dreal( 1d0 / r )

!    x = dreal( 1d0 / r )
!    y = dimag( 1d0 / r )
!    phi = datan2(y, x)
!    w_r_recp_abs = (x*x+y*y)**0.25 * cdexp(eye*0.5*(phi+k*PI))
!    w_r_recp_abs = (x*x+y*y)**0.25 * (dcos(0.5*(phi+k*PI)) + eye*dsin(0.5*(phi)))

    return
end function sub_int_r_reciprocal_abs
