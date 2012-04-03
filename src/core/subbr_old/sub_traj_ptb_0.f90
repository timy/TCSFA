#include '../../include/inc_field.h'

double complex function sub_traj_x_0( t, ts )
    use mod_p0

    implicit none;
    double complex, intent(in):: t, ts
    double complex, external:: PULSE_ALPHA_X

    sub_traj_x_0 = PULSE_ALPHA_X( t ) + p0_x * t - &
          dreal( PULSE_ALPHA_X( ts ) + p0_x * ts );

    return;
end function sub_traj_x_0


double complex function sub_traj_z_0( t, ts )
    use mod_p0

    implicit none;
    double complex, intent(in):: t, ts;
    double complex, external:: PULSE_ALPHA_Z;

    sub_traj_z_0 = PULSE_ALPHA_Z( t ) + p0_z * t - &
          dreal( PULSE_ALPHA_Z( ts ) + p0_z * ts );

    return;
end function sub_traj_z_0


double complex function sub_traj_vx_0( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double complex, external:: PULSE_A_X;

    sub_traj_vx_0 = p0_x + PULSE_A_X( t );

    return;
end function sub_traj_vx_0


double complex function sub_traj_vz_0( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double complex, external:: PULSE_A_Z;

    sub_traj_vz_0 = p0_z + PULSE_A_Z( t );

    return;
end function sub_traj_vz_0
