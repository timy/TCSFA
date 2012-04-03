#include '../../include/inc_atom.h'

! ////////////////////////////////////////////////////////////////////////////////
! the action correction due to the -1/r term of  the sub-barrier trajectory 
double complex function action_W_sub_r_rcpr( ts )
    implicit none
    double complex:: ts;
    double precision:: ti
    integer, parameter:: n_step = 2000
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: simpson_sub, integrand_sub_r_rcpr_1st
 
    ti = dimag(ts)    
    call generate_perturb_vx(ts)
    call generate_perturb_x(ts)
    action_W_sub_r_rcpr = eye * simpson_sub(0d0, ti, n_step, ts, integrand_sub_r_rcpr_1st)
    return
end function action_W_sub_r_rcpr

! ////////////////////////////////////////////////////////////////////////////////
! the action correction based on Sergey's perturbation method
double complex function action_W_sub_sergey( ts )
    implicit none
    double complex:: ts;
    double precision:: n_star
    double precision:: ti
    integer, parameter:: n_step = 2000
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: simpson_sub, integrand_sub_r_rcpr_0th
    double precision, parameter:: epsilon = 0.0001d0

    ti = dimag(ts)
    ! the 1st term (sergey)
    n_star = ATOM_CHARGE_Z / dsqrt(2d0*IONIZATION_IP)
    action_W_sub_sergey = -eye * n_star * dlog(2d0*IONIZATION_IP*ti)
    ! plus the 2nd term
    action_W_sub_sergey = action_W_sub_sergey + &
          eye*simpson_sub( 0d0, ti-epsilon, n_step, ts, integrand_sub_r_rcpr_0th )
    return

end function action_W_sub_sergey

! ////////////////////////////////////////////////////////////////////////////////
! the integrand
double complex function integrand_sub_r_rcpr_1st( t, ts ) result(integrand)
    implicit none
    double precision:: t, t0, ti
    double complex:: ts, z, x, r, tau
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: sub_traj_z_0, sub_traj_x_0, sub_traj_x_1

    t0 = dreal(ts)
    ti = dimag(ts)
    tau = t0 + eye*t
    z = sub_traj_z_0(tau, ts)
    x = sub_traj_x_0(tau, ts) + sub_traj_x_1(t)
    r = cdsqrt(z*z+x*x)
    integrand = - ATOM_CHARGE_Z / r
    return
end function integrand_sub_r_rcpr_1st

! ////////////////////////////////////////////////////////////////////////////////
! the integrand for direct 0th-order trajectory
double complex function integrand_sub_r_rcpr_0th( t, ts ) result(integrand)
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double precision:: t0, ti
    double complex:: z, x, r, tau
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: sub_traj_z_0, sub_traj_x_0

    t0 = dreal(ts)
    ti = dimag(ts)
    tau = t0 + eye*t
    z = sub_traj_z_0(tau, ts)
    x = sub_traj_x_0(tau, ts)
    r = cdsqrt(z*z+x*x)
    integrand = - ATOM_CHARGE_Z / r
    return
end function integrand_sub_r_rcpr_0th
