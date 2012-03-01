#include '../../include/inc_atom.h'

! ////////////////////////////////////////////////////////////////////////////////
! Two modules for the z-direction calculation
module mod_interp_vz
    implicit none
    integer, parameter:: n_samp = 40
    double precision:: samp_dt
    double precision:: samp_t(n_samp)
    double complex:: samp_vz(n_samp), b_vz(n_samp), c_vz(n_samp), d_vz(n_samp)
end module mod_interp_vz

module mod_interp_z
    implicit none
    integer, parameter:: n_samp = 40
    double precision:: samp_dt
    double precision:: samp_t(n_samp)
    double complex:: samp_z(n_samp), b_z(n_samp), c_z(n_samp), d_z(n_samp)
end module mod_interp_z

! Two modules for the x-direction calculation
module mod_interp_vx
    implicit none
    integer, parameter:: n_samp = 40
    double precision:: samp_dt
    double precision:: samp_t(n_samp)
    double complex:: samp_vx(n_samp), b_vx(n_samp), c_vx(n_samp), d_vx(n_samp)
end module mod_interp_vx

module mod_interp_x
    implicit none
    integer, parameter:: n_samp = 40
    double precision:: samp_dt
    double precision:: samp_t(n_samp)
    double complex:: samp_x(n_samp), b_x(n_samp), c_x(n_samp), d_x(n_samp)
end module mod_interp_x

! ////////////////////////////////////////////////////////////////////////////////
! the force in the z-direction with perturbation method
double complex function perturb_fz(t, ts)
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex, external:: sub_traj_x_0, sub_traj_z_0
    double complex:: x_0, z_0
    double precision:: t0
    t0 = dreal(ts)
    x_0 = sub_traj_x_0( dcmplx(t0, t), ts )
    z_0 = sub_traj_z_0( dcmplx(t0, t), ts )
    perturb_fz = - ATOM_CHARGE_Z * z_0 / ( x_0 * x_0  + z_0 * z_0 ) ** (1.5)
    return
end function perturb_fz

! the force in the x-direction with perturbation method
double complex function perturb_fx(t, ts)
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex, external:: sub_traj_x_0, sub_traj_z_0
    double complex:: x_0, z_0
    double precision:: t0
    t0 = dreal(ts)
    x_0 = sub_traj_x_0( dcmplx(t0, t), ts )
    z_0 = sub_traj_z_0( dcmplx(t0, t), ts )
    perturb_fx = - ATOM_CHARGE_Z * x_0 / ( x_0 * x_0  + z_0 * z_0 ) ** (1.5)
    return
end function perturb_fx

! ////////////////////////////////////////////////////////////////////////////////
! generate the interpolation info for the sub_traj_vz with perturbation method
subroutine generate_perturb_vz(ts)
    use mod_interp_vz, only: n_samp, samp_dt, samp_t, samp_vz, &
          b_vz, c_vz, d_vz
    implicit none
    double complex, intent(in):: ts
    double precision:: ti
    double complex:: c1
    double complex, external:: perturb_fz, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i
    ti = dimag(ts)
    samp_dt = (ti - 0d0) / (n_samp - 1d0)
    c1 = - eye * dimag( eye * simpson_sub( ti, 0d0, 1000, ts, perturb_fz ) )

    forall(i = 1:n_samp) samp_t(i) = 0d0 + samp_dt * (i - 1d0)
    do i = 1, n_samp
        samp_vz(i) = eye * simpson_sub(ti, samp_t(i), 1000, ts, perturb_fz) + c1
    end do
    call spline(n_samp, samp_t, samp_vz, b_vz, c_vz, d_vz)
    return
end subroutine generate_perturb_vz

! generate the interpolation info for the sub_traj_vx with perturbation method
subroutine generate_perturb_vx(ts)
    use mod_interp_vx, only: n_samp, samp_dt, samp_t, samp_vx, &
          b_vx, c_vx, d_vx
    implicit none
    double complex, intent(in):: ts
    double precision:: ti
    double complex:: c1
    double complex, external:: perturb_fx, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i
    ti = dimag(ts)
    samp_dt = (ti - 0d0) / (n_samp - 1d0)
    c1 = - eye * dimag( eye * simpson_sub( ti, 0d0, 1000, ts, perturb_fx ) )

    forall(i = 1:n_samp) samp_t(i) = 0d0 + samp_dt * (i - 1d0)
    do i = 1, n_samp
        samp_vx(i) = eye * simpson_sub(ti, samp_t(i), 1000, ts, perturb_fx) + c1
    end do
    call spline(n_samp, samp_t, samp_vx, b_vx, c_vx, d_vx)
    return
end subroutine generate_perturb_vx

! ////////////////////////////////////////////////////////////////////////////////
! the final interpolation function for the sub_traj_vz with perturbation method
double complex function sub_traj_vz_1(t, ts)
    use mod_interp_vz, only: samp_dt, samp_t, samp_vz, b_vz, c_vz, d_vz
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    integer:: index
    index = ceiling( t / samp_dt )
    if(index .eq. 0) index = 1
    sub_traj_vz_1 = samp_vz(index) + b_vz(index) * ( t - samp_t(index) ) + &
         c_vz(index) * ( t - samp_t(index) )**2 + d_vz(index) * ( t - samp_t(index) )**3
end function sub_traj_vz_1

! the final interpolation function for the sub_traj_vx with perturbation method
double complex function sub_traj_vx_1(t, ts)
    use mod_interp_vx, only: samp_dt, samp_t, samp_vx, b_vx, c_vx, d_vx
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    integer:: index
    index = ceiling( t / samp_dt )
    if(index .eq. 0) index = 1
    sub_traj_vx_1 = samp_vx(index) + b_vx(index) * ( t - samp_t(index) ) + &
          c_vx(index) * ( t - samp_t(index) )**2 + d_vx(index) * ( t - samp_t(index) )**3
end function sub_traj_vx_1

! ////////////////////////////////////////////////////////////////////////////////
! generate the interpolation info for the sub_traj_z with perturbation method
subroutine generate_perturb_z(ts)
    use mod_interp_z, only: n_samp, samp_dt, samp_t, samp_z, b_z, c_z, d_z
    implicit none
    double complex, intent(in):: ts
    double precision:: ti
    double complex:: c2
    double complex, external:: sub_traj_vz_1, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i
    ti = dimag(ts)
    samp_dt = (ti - 0d0) / (n_samp - 1d0)
    c2 = - eye*dimag( eye * simpson_sub( ti, 0d0, 1000, ts, sub_traj_vz_1 ) )
    forall(i = 1:n_samp) samp_t(i) = 0d0 + samp_dt * (i - 1d0)
    do i = 1, n_samp
        samp_z(i) = eye * simpson_sub(ti, samp_t(i), 1000, ts, sub_traj_vz_1) + c2
    end do
    call spline(n_samp, samp_t, samp_z, b_z, c_z, d_z)
    return
end subroutine generate_perturb_z

! generate the interpolation info for the sub_traj_x with perturbation method
subroutine generate_perturb_x(ts)
    use mod_interp_x, only: n_samp, samp_dt, samp_t, samp_x, b_x, c_x, d_x
    implicit none
    double complex, intent(in):: ts
    double precision:: ti
    double complex:: c2
    double complex, external:: sub_traj_vx_1, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i
    ti = dimag(ts)
    samp_dt = (ti - 0d0) / (n_samp - 1d0)
    c2 = - eye*dimag( eye * simpson_sub( ti, 0d0, 1000, ts, sub_traj_vx_1 ) )
    forall(i = 1:n_samp) samp_t(i) = 0d0 + samp_dt * (i - 1d0)
    do i = 1, n_samp
        samp_x(i) = eye * simpson_sub(ti, samp_t(i), 1000, ts, sub_traj_vx_1) + c2
    end do
    call spline(n_samp, samp_t, samp_x, b_x, c_x, d_x)
    return
end subroutine generate_perturb_x
! ////////////////////////////////////////////////////////////////////////////////
! the final interpolation function for the sub_traj_z with perturbation method
double complex function sub_traj_z_1(t)
    use mod_interp_z, only: samp_dt, samp_t, samp_z, b_z, c_z, d_z
    implicit none
    double precision, intent(in):: t
    integer:: index
    index = ceiling( t / samp_dt )
    if(index .eq. 0) index = 1
    sub_traj_z_1 = samp_z(index) + b_z(index) * ( t - samp_t(index) ) + &
         c_z(index) * ( t - samp_t(index) )**2 + d_z(index) * ( t - samp_t(index) )**3
end function sub_traj_z_1

! the final interpolation function for the sub_traj_x with perturbation method
double complex function sub_traj_x_1(t)
    use mod_interp_x, only: samp_dt, samp_t, samp_x, b_x, c_x, d_x
    implicit none
    double precision, intent(in):: t
    integer:: index
    index = ceiling( t / samp_dt )
    if(index .eq. 0) index = 1
    sub_traj_x_1 = samp_x(index) + b_x(index) * ( t - samp_t(index) ) + &
         c_x(index) * ( t - samp_t(index) )**2 + d_x(index) * ( t - samp_t(index) )**3
end function sub_traj_x_1
! ////////////////////////////////////////////////////////////////////////////////
! plot the sub-barrier trajectories with the perturbation method
subroutine plot_sub_traj_ptb_1_interp(ts)
    implicit none
    double complex, intent(in):: ts
    integer:: i
    integer, parameter:: n_samp = 200
    double precision:: ti, dt, t
    double complex, external:: perturb_fz, sub_traj_vz_1, sub_traj_z_1
    double complex, external:: perturb_fx, sub_traj_vx_1, sub_traj_x_1

    ti = dimag(ts)-1e-6
    dt = ( ti - 0d0 ) / ( n_samp - 1d0 )

    open(101, file="dat/sub_traj_interp.dat")
    do i = 1, n_samp
        t = 0d0 + ( i - 1d0 ) * dt
        write(101, '(13(e15.8,1x))'), t, perturb_fz(t, ts), perturb_fx(t, ts), &
             sub_traj_vz_1(t, ts), sub_traj_vx_1(t, ts), &
             sub_traj_z_1(t), sub_traj_x_1(t)
    end do
    close(101)

    return
end subroutine plot_sub_traj_ptb_1_interp

! ////////////////////////////////////////////////////////////////////////////////
! the test module for generate_perturb_vx and generate_perturb_x
! plot the sample points used for interpolation

subroutine plot_sub_traj_ptb_1_samp(ts)
    use mod_interp_vx, only: n_samp
    implicit none
    double complex, intent(in):: ts
    double precision:: ti, dt, t(n_samp)
    double complex:: cz1, cz2, vz(n_samp), z(n_samp)
    double complex:: cx1, cx2, vx(n_samp), x(n_samp)
    double complex, external:: perturb_fz, sub_traj_vz_1, simpson_sub
    double complex, external:: perturb_fx, sub_traj_vx_1
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i

    ti = dimag(ts)
    dt = (ti - 0d0) / (n_samp - 1d0)
    cz1 = - eye * dimag( eye * simpson_sub( ti, 0d0, 1000, ts, perturb_fz ) )
    cz2 = - eye * dimag( simpson_sub( ti, 0d0, 1000, ts, sub_traj_vz_1 ) )
    cx1 = - eye * dimag( eye * simpson_sub( ti, 0d0, 1000, ts, perturb_fx ) )
    cx2 = - eye * dimag( simpson_sub( ti, 0d0, 1000, ts, sub_traj_vx_1 ) )
    forall(i = 1:n_samp) t(i) = 0d0 + dt * (i - 1d0)
    do i = 1, n_samp
        vz(i) = eye * simpson_sub(ti, t(i), 1000, ts, perturb_fz) + cz1
        z(i) = simpson_sub(ti, t(i), 1000, ts, sub_traj_vz_1) + cz2
        vx(i) = eye * simpson_sub(ti, t(i), 1000, ts, perturb_fx) + cx1
        x(i) = simpson_sub(ti, t(i), 1000, ts, sub_traj_vx_1) + cx2
    end do

    open(101, file="dat/sub_traj_samp.dat")
    do i = 1, n_samp
        write(101, '(9(e15.8,1x))'), t(i), vz(i), z(i), vx(i), x(i)
    end do
    close(101)

    return
end subroutine plot_sub_traj_ptb_1_samp

subroutine plot_sub_traj_ptb_1(ts)
    implicit none
    double complex, intent(in):: ts
    call plot_sub_traj_ptb_1_interp(ts)
    call plot_sub_traj_ptb_1_samp(ts)
    stop
end subroutine plot_sub_traj_ptb_1
