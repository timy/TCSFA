#include '../include/inc_atom.h'
#include '../include/inc_misc.h'
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
    double complex, external:: pulse_A_x, pulse_A_z;

    
    SPE = ( p0_x + pulse_A_x(t) )**2 + ( p0_z + pulse_A_z(t) )**2 + 2d0 * Ip;

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
    double complex, external:: pulse_alpha_x

    im_traj_x = pulse_alpha_x( t ) + p0_x * t - &
          dreal( pulse_alpha_x( ts ) + p0_x * ts );

    return;
end function im_traj_x


double complex function im_traj_z( t, ts )
    use mod_p0

    implicit none;
    double complex, intent(in):: t, ts;
    double complex, external:: pulse_alpha_z;

    im_traj_z = pulse_alpha_z( t ) + p0_z * t - &
          dreal( pulse_alpha_z( ts ) + p0_z * ts );

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


! double complex function action_W_im( ts )
!     use mod_p0
!     use mod_pulse, only: om, E0, nc, xi, Pi
    
!     implicit none;
!     double complex, intent(in):: ts;
!     double complex, external:: v2_integrand, cedint;
!     double complex, external:: pulse_alpha_z, pulse_alpha2_z
!     double complex:: t0;
!     double precision, parameter:: Ip = IONIZATION_IP
!     double precision:: p0x, p0z
! #ifdef IM_PLOT_INTEGRAND
!     integer:: i
!     double complex:: t
! #endif

!     t0 = dcmplx( dreal(ts), 0d0 );
!     p0x = p0_x;
!     p0z = p0_z;
!     if ( dimag(ts)+1d-6 <= 0d0 ) stop 'error: Im[t] < 0';

! #ifdef IM_PLOT_INTEGRAND
!      open(IM_PLOT_FILE_ID, file=IM_PLOT_FILE_NAME)
!      do i = 1, IM_PLOT_N_PTS
!          t = dcmplx( dreal(ts), ( dimag(ts) + IM_PLOT_OFFSET ) / (IM_PLOT_N_PTS-1) * (i - 1) );
!          write(IM_PLOT_FILE_ID, '(3(e15.8,2x))'), dimag(t), v2_integrand( t ) + Ip;
!      end do
!      close(IM_PLOT_FILE_ID)
! #endif ! IM_PLOT_INTEGRAND
     
!     action_w_im = (0.5d0 * (p0x*p0x+p0z*p0z) + Ip ) * ( t0 - ts ) &
!           + p0z * ( pulse_alpha_z(t0) - pulse_alpha_z(ts) ) &
!           + 0.5d0 * ( pulse_alpha2_z(t0) - pulse_alpha2_z(ts) )

!     return;
! end function action_W_im


