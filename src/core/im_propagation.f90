#include '../include/inc_atom.h'
#include '../include/inc_misc.h'
#include '../include/inc_field.h'
#ifdef MISC_PLOT
#include '../include/inc_plot_im_integrand.h'
#endif

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

double complex function solve_ts_from_p0( ts0, ierr ) result(ts)
    use mod_p0;

    implicit none;

    double complex, intent(in):: ts0;
    double complex, external:: SPE
    double complex, external:: find_root
    integer:: ierr

#if MISC_PRINT > 3
    double complex, external:: im_traj_x, im_traj_z, im_traj_vx, im_traj_vz
#endif

    ts = find_root( ts0, SPE, CRF_STEP_LENGTH, ierr );
    if( ierr > 0 ) then
        print*, 'root not found!'
        return;
    end if

#if MISC_PRINT > 3
    print*, 'ts:', ts
    print*, 'x(ts):', im_traj_x( ts, ts );
    print*, 'z(ts):', im_traj_z( ts, ts );
    print*, 'vx(ts):', im_traj_vx( ts );
    print*, 'vz(ts):', im_traj_vz( ts );
#endif

end function solve_ts_from_p0




double complex function SPE( t ) ! saddle point equation
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double precision, parameter:: Ip = IONIZATION_IP;
    double complex, external:: PULSE_A_X, PULSE_A_Z;

    
    SPE = ( p0_x + PULSE_A_X(t) )**2 + ( p0_z + PULSE_A_Z(t) )**2 + 2d0 * Ip;

    return;
end function SPE





! calculate action W for imaginary time part with numerical integration
double complex function action_W_im_num( ts )
    
    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: v2_integrand, cedint;
    double complex:: t0;
    double precision, parameter:: Ip = IONIZATION_IP
#ifdef IM_PLOT_INTEGRAND
    integer:: i
    double complex:: t
#endif

    t0 = dcmplx( dreal(ts), 0d0 );
    if ( dimag(ts)+1d-6 <= 0d0 ) stop 'error: Im[t] < 0';

#ifdef IM_PLOT_INTEGRAND
     open(IM_PLOT_FILE_ID, file=IM_PLOT_FILE_NAME)
     do i = 1, IM_PLOT_N_PTS
         t = dcmplx( dreal(ts), ( dimag(ts) + IM_PLOT_OFFSET ) / (IM_PLOT_N_PTS-1) * (i - 1) );
         write(IM_PLOT_FILE_ID, '(3(e15.8,2x))'), dimag(t), v2_integrand( t ) + Ip;
     end do
     close(IM_PLOT_FILE_ID)
#endif ! IM_PLOT_INTEGRAND

    action_W_im_num = cedint( v2_integrand, 1000, ts, t0 );
    action_W_im_num = action_W_im_num + ( t0 - ts ) * Ip
    
    return;
end function action_W_im_num



! second-order derivative of the action W
double complex function action_DDW( ts )
    
    implicit none
    double complex, intent(in):: ts
    double complex:: vzs, vxs, azs, axs
    double complex, external:: PULSE_E_Z, PULSE_E_X, im_traj_vz, im_traj_vx

    vzs = im_traj_vz( ts );
    vxs = im_traj_vx( ts );
    azs = - PULSE_E_Z( ts );
    axs = - PULSE_E_X( ts );
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
    double complex, external:: PULSE_ALPHA_X

    im_traj_x = PULSE_ALPHA_X( t ) + p0_x * t - &
          dreal( PULSE_ALPHA_X( ts ) + p0_x * ts );

    return;
end function im_traj_x


double complex function im_traj_z( t, ts )
    use mod_p0

    implicit none;
    double complex, intent(in):: t, ts;
    double complex, external:: PULSE_ALPHA_Z;

    im_traj_z = PULSE_ALPHA_Z( t ) + p0_z * t - &
          dreal( PULSE_ALPHA_Z( ts ) + p0_z * ts );

    return;
end function im_traj_z


double complex function im_traj_vx( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double complex, external:: PULSE_A_X;

    im_traj_vx = p0_x + PULSE_A_X( t );

    return;
end function im_traj_vx


double complex function im_traj_vz( t )
    use mod_p0;

    implicit none;
    double complex, intent(in):: t;
    double complex, external:: PULSE_A_Z;

    im_traj_vz = p0_z + PULSE_A_Z( t );

    return;
end function im_traj_vz
