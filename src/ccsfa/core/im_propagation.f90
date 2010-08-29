#include '../include/inc_atom.h'

!#define IM_PROP_PRNT_DATA

#define CRF_STEP_LENGTH 0.1d0


! in this file, saddle_point equation should be defined


module mod_p0
    implicit none;
    double precision:: p0_x;
    double precision:: p0_z;
end module mod_p0


subroutine set_p0( p0_x_, p0_z_ )
    use mod_p0
    
    implicit none
    double precision, intent(in):: p0_x_, p0_z_

    p0_x = p0_x_;
    p0_z = p0_z_;
    
    return;
end subroutine set_p0



double complex function solve_ts_from_p0( ts0 ) result(ts)
    use mod_p0;

    implicit none;

    double complex, intent(in):: ts0;
    double complex, external:: SPE
    double complex:: z0
    double precision:: hs
    integer:: flag
    double complex, external:: find_root
    integer:: ierr

#ifdef IM_PROP_PRNT_DATA
    double complex, external:: im_traj_x, im_traj_z, im_traj_vx, im_traj_vz
#endif

    flag = 0; ! used for further extension

    ts = find_root( ts0, SPE, CRF_STEP_LENGTH, flag, ierr );
    if( ierr > 0 ) then
        print*, 'root not found!'
        return;
    end if

#ifdef IM_PROP_PRNT_DATA
    print*, 'ts:', ts
    print*, 'x(ts):', im_traj_x( ts, ts );
    print*, 'z(ts):', im_traj_z( ts, ts );
    print*, 'vx(ts):', im_traj_vx( ts );
    print*, 'vz(ts):', im_traj_vz( ts );
#endif

end function solve_ts_from_p0




double complex function SPE( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double precision, parameter:: Ip = IONIZATION_IP;
    double complex, external:: pulse_A_x, pulse_A_z;

    
    SPE = ( p0_x + pulse_A_x(t) )**2 + ( p0_z + pulse_A_z(t) )**2 + 2d0 * Ip;

    return;
end function SPE






double complex function action_W_im( ts )
    
    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: v2_integrand, cedint;
    double complex:: t0;
    double precision, parameter:: Ip = IONIZATION_IP

    t0 = dcmplx( dreal(ts), 0d0 );
    if ( dimag(ts)+1d-6 <= 0d0 ) stop 'error: Im[t] < 0';

    action_W_im = cedint( v2_integrand, 1000, ts, t0 );
    action_W_im = action_W_im + ( t0 - ts ) * Ip
    
    return;
end function action_W_im




double complex function action_DDW( ts )
    
    implicit none
    double complex, intent(in):: ts
    double complex:: vzs, vxs, azs, axs
    double complex, external:: pulse_E_z, pulse_E_x, im_traj_vz, im_traj_vx

    vzs = im_traj_vz( ts );
    vxs = im_traj_vx( ts );
    azs = - pulse_E_z( ts );
    axs = - pulse_E_x( ts );
    action_DDW = vzs * azs + vxs * axs;
    
    return;
end function action_DDW


double complex function v2_integrand( t )

    implicit none
    double complex, intent(in):: t
    double complex, external:: im_traj_x, im_traj_z, im_traj_vx, im_traj_vz
    double complex:: vz, vx, v2
    
    vz = im_traj_vz( t );
    vx = im_traj_vx( t );
    v2 = vz * vz + vx * vx;

    v2_integrand = 0.5d0 * v2
    
    return;
end function v2_integrand


















double complex function im_traj_x( t, ts )
    use mod_p0

    implicit none;
    double complex, intent(in):: t, ts
    double complex, external:: alpha_x

    im_traj_x = alpha_x( t ) + p0_x * t - &
          dreal( alpha_x( ts ) + p0_x * ts );

    return;
end function im_traj_x


double complex function im_traj_z( t, ts )
    use mod_p0

    implicit none;
    double complex, intent(in):: t, ts;
    double complex, external:: alpha_z;

    im_traj_z = alpha_z( t ) + p0_z * t - &
          dreal( alpha_z( ts ) + p0_z * ts );

    return;
end function im_traj_z


double complex function im_traj_vx( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double complex, external:: pulse_A_x;

    im_traj_vx = p0_x + pulse_A_x( t );

    return;
end function im_traj_vx


double complex function im_traj_vz( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double complex, external:: pulse_A_z;

    im_traj_vz = p0_z + pulse_A_z( t );

    return;
end function im_traj_vz
