#include '../../include/inc_atom.h'

! ////////////////////////////////////////////////////////////////////////////////
! the action correction due to the Coulomb influence on the sub-barrier trajectory 
double complex function action_S_sub( ts )
    implicit none
    double complex:: ts;
    double precision:: ti
    integer, parameter:: n_step = 2000
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: simpson_sub, integrand_sub
 
    ti = dimag(ts)    
    call generate_perturb_vx(ts)
    call generate_perturb_x(ts)
    action_S_sub = eye * simpson_sub(0d0, ti, n_step, ts, integrand_sub)
!    print*, "action_S_sub with perturbation", action_S_sub
    
    return

end function action_S_sub

! ////////////////////////////////////////////////////////////////////////////////
! the action correction based on Sergey's perturbation method
double complex function action_S_sub_sergey( ts )
    implicit none
    double complex:: ts;
    double precision:: n_star
    double precision:: ti
    integer, parameter:: n_step = 2000
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: simpson_sub, integrand_sub
    double precision, parameter:: epsilon = 0.0001d0

    ti = dimag(ts)
    ! the 1st term (sergey)
    n_star = ATOM_CHARGE_Z / dsqrt(2d0*IONIZATION_IP)
    action_S_sub_sergey = -eye * n_star * dlog(2d0*IONIZATION_IP*ti)
    ! plus the 2nd term
    action_S_sub_sergey = action_S_sub_sergey + &
          eye*simpson_sub( 0d0, ti-epsilon, n_step, ts, integrand_sub )

! obviously, simpson_sub performs better than integration_sub
!       print*, "trapisod:", trapezoid_sub( 0d0, ti-epsilon, n_step, ts )
!       print*, "simpson:", simpson_sub( 0d0, ti-epsilon, n_step, ts )
!       print*, "action_S_sub:", action_S_sub
    
    return

end function action_S_sub_sergey

! ////////////////////////////////////////////////////////////////////////////////
! the integrand
double complex function integrand_sub( t, ts )
!    use mod_p0
    implicit none
    double precision:: t, t0, ti
    double complex:: ts, z, x, r, tau
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: im_traj_z, im_traj_x, interp_perturb_x

    t0 = dreal(ts)
    ti = dimag(ts)
    tau = t0 + eye*t

! my coulomb correction with 1st order perturbation expansion in the x-direction
    z = im_traj_z(tau, ts)
    x = interp_perturb_x(t)
    r = cdsqrt(z*z+x*x)
    integrand_sub = - ATOM_CHARGE_Z / r
! end my coulomb correction
    return
end function integrand_sub
