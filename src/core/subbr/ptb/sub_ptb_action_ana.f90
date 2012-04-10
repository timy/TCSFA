#include '../../../include/inc_atom.h'
#include '../../../include/inc_field.h'

! the sub-barrier action for the 0th-order trajectory without -1/r term
! for the linearly polarized field along the z-axis
double complex function action_W_sub_0_0_0_ana( ts ) result(w_sub)
    use mod_p0, only: p0_x, p0_z
    implicit none
    double complex, intent(in):: ts
    double complex, external:: PULSE_ALPHA_Z, PULSE_ALPHA2_Z
    double complex:: t0

    t0 = dcmplx( dreal(ts), 0d0 )
    w_sub = &
         (0.5d0 * (p0_z*p0_z+p0_x*p0_x) + IONIZATION_IP ) * ( ts - t0 ) &
         + p0_z * ( PULSE_ALPHA_Z(ts) - PULSE_ALPHA_Z(t0) ) &
         + 0.5d0 * ( PULSE_ALPHA2_Z(ts) - PULSE_ALPHA2_Z(t0) )
    return;
end function action_W_sub_0_0_0_ana

! the sub-barrier action for the 0th-order trajectory with -1/r term
! for the linearly polarized field along the z-axis
double complex function action_w_sub_1_0_0_ana(ts) result(w_sub)
    use mod_p0, only: p0_x, p0_z
    implicit none
    double complex, intent(in):: ts
    double complex, external:: PULSE_ALPHA_Z, PULSE_ALPHA2_Z
    double complex, external:: sub_int_r_reciprocal, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: t0, w_r_recp
    double precision:: ti

    t0 = dcmplx( dreal(ts), 0d0 )
    ti = dimag(ts)
    w_sub = &
         (0.5d0 * (p0_z*p0_z+p0_x*p0_x) + IONIZATION_IP ) * ( ts - t0 ) &
         + p0_z * ( PULSE_ALPHA_Z(ts) - PULSE_ALPHA_Z(t0) ) &
         + 0.5d0 * ( PULSE_ALPHA2_Z(ts) - PULSE_ALPHA2_Z(t0) )

    w_r_recp = - eye * simpson_sub(ti, 0d0, 200, ts, sub_int_r_reciprocal)
    w_sub = w_sub + w_r_recp

    return
end function action_w_sub_1_0_0_ana

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
    return
end function sub_int_r_reciprocal
