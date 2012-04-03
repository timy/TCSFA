#include '../../include/inc_field.h'
#include '../../include/inc_atom.h'

! ////////////////////////////////////////////////////////////////////////////////
! the sub-barrier action for the 0th-order trajectory without -1/r term
! for the linearly polarized field along the z-axis
double complex function action_W_sub_0_traj_0( ts )
    use mod_p0, only: p0_x, p0_z
    implicit none
    double complex, intent(in):: ts
    double complex, external:: PULSE_ALPHA_Z, PULSE_ALPHA2_Z
    double complex:: t0

    t0 = dcmplx( dreal(ts), 0d0 )
    action_W_sub_0_traj_0 = &
         (0.5d0 * (p0_z*p0_z+p0_x*p0_x) + IONIZATION_IP ) * ( ts - t0 ) &
         + p0_z * ( PULSE_ALPHA_Z(ts) - PULSE_ALPHA_Z(t0) ) &
         + 0.5d0 * ( PULSE_ALPHA2_Z(ts) - PULSE_ALPHA2_Z(t0) )
    return;
end function action_W_sub_0_traj_0
