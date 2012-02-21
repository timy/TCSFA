#include '../../include/inc_atom.h'

! ////////////////////////////////////////////////////////////////////////////////
! the general form of action for sub-barrier regime
double complex function action_W_sub( ts ) result(action)
    implicit none
    double complex:: ts;
    double precision:: ti
    integer, parameter:: n_step = 2000
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: simpson_sub, integrand_sub_1st
 
    ti = dimag(ts)    
    call generate_perturb_vx(ts)
    call generate_perturb_x(ts)
    action = eye * simpson_sub(0d0, ti, n_step, ts, integrand_sub_1st)
    return
end function action_W_sub

! ////////////////////////////////////////////////////////////////////////////////
! the integrand
double complex function integrand_sub_1st( t, ts ) result(integrand)
    implicit none
    double precision:: t, t0, ti
    double complex:: ts, z, x, r, tau, vz, vx, v2
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: sub_traj_z_0, sub_traj_x_0, sub_traj_x_1, &
          sub_traj_vz_0, sub_traj_vx_0, sub_traj_vx_1

    t0 = dreal(ts)
    ti = dimag(ts)
    tau = t0 + eye*t
    vz = sub_traj_vz_0( tau )
    vx = sub_traj_vx_0( tau ) !+ sub_traj_vx_1( t, ts )
    v2 = vz * vz + vx * vx
    z = sub_traj_z_0(tau, ts)
    x = sub_traj_x_0(tau, ts) + sub_traj_x_1(t)
    r = cdsqrt(z*z+x*x)
    integrand = 0.5 * v2 - ATOM_CHARGE_Z / r + IONIZATION_IP
    return
end function integrand_sub_1st
