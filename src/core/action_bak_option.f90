double complex function integrand_sub( t, ts )
    use mod_p0
    implicit none
    double precision:: t, t0, ti, n_star, p
    double complex:: ts, z, x, r, tau
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex, external:: im_traj_z, im_traj_x, interp_perturb_x
    double complex, external:: pulse_alpha_z, pulse_alpha_x

    print*, "integrand_sub"

    n_star = ATOM_CHARGE_Z / dsqrt(2d0*IONIZATION_IP)
    t0 = dreal(ts)
    ti = dimag(ts)
    tau = t0 + eye*t
!    z = im_traj_z( tau, ts )
!    x = im_traj_x( tau, ts )

!    r = cdsqrt(z*z+x*x)
!    r = cdsqrt(z*z)

! new strategy 1 from dieter
!    z = im_traj_z( tau, ts )
!    r = cdsqrt( z*z )
!new strategy end

! new strategy 2  from dieter
!    p = dsqrt(p0_z*p0_z+p0_x*p0_x)
!    z = pulse_alpha_z(tau) + p*tau - dreal(pulse_alpha_z(ts) + p*ts)
!    r = cdsqrt(z*z)
! new strategy end

! ! sergey's original correction
!     z = pulse_alpha_z(tau) - pulse_alpha_z(ts) + p0_z*(tau-ts)
!     x = p0_x*(tau-ts)
!     r = cdsqrt(z*z+x*x)
!     integrand_sub = n_star / (ti - t) - ATOM_CHARGE_Z / r
! ! end sergey's original correction

! my coulomb correction
!     z = im_traj_z(tau, ts)
!     x = im_traj_x(tau, ts)
!     r = cdsqrt(z*z+x*x)
!     integrand_sub = - ATOM_CHARGE_Z / r
! end my coulomb correcttion

! my coulomb correction with 0th order perturbation expansion in the x-direction
    z = im_traj_z(tau, ts)
    x = interp_perturb_x(t)
    r = cdsqrt(z*z+x*x)
    integrand_sub = - ATOM_CHARGE_Z / r
! end my coulomb correction
    return
end function integrand_sub
