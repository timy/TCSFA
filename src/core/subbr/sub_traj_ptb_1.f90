#include '../../include/inc_atom.h'

! ////////////////////////////////////////////////////////////////////////////////
! Two modules for the following calculation
module mod_interp_vx
    implicit none
    integer, parameter:: n_samp = 40
    double precision:: samp_dt
    double precision:: samp_t(n_samp)
    double complex:: samp_vx(n_samp), b(n_samp), c(n_samp), d(n_samp)
end module mod_interp_vx

module mod_interp_x
    implicit none
    integer, parameter:: n_samp = 40
    double precision:: samp_dt
    double precision:: samp_t(n_samp)
    double complex:: samp_x(n_samp), b(n_samp), c(n_samp), d(n_samp)
end module mod_interp_x

! ////////////////////////////////////////////////////////////////////////////////
! the force in the x-direction with perturbation method
double complex function perturb_fx(t, ts)
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex, external:: im_traj_x, im_traj_z
    double complex:: x_0, z_0
    double precision:: t0
    t0 = dreal(ts)
    x_0 = im_traj_x( dcmplx(t0, t), ts )
    z_0 = im_traj_z( dcmplx(t0, t), ts )
    perturb_fx = - ATOM_CHARGE_Z * x_0 / ( x_0 * x_0  + z_0 * z_0 ) ** (1.5)
    return
end function perturb_fx

! ////////////////////////////////////////////////////////////////////////////////
! generate the interpolation info for the im_traj_vx with perturbation method
subroutine generate_perturb_vx(ts)
    use mod_interp_vx, only: n_samp, samp_dt, samp_t, samp_vx, b, c, d
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
    call spline(n_samp, samp_t, samp_vx, b, c, d)
    return
end subroutine generate_perturb_vx

! ////////////////////////////////////////////////////////////////////////////////
! the final interpolation function for the im_traj_vx with perturbation method
double complex function interp_perturb_vx(t, ts)
    use mod_interp_vx, only: samp_dt, samp_t, samp_vx, b, c, d
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    integer:: index
    index = ceiling( t / samp_dt )
    if(index .eq. 0) index = 1
    interp_perturb_vx = samp_vx(index) + b(index) * ( t - samp_t(index) ) + &
         c(index) * ( t - samp_t(index) )**2 + d(index) * ( t - samp_t(index) )**3
end function interp_perturb_vx

! ////////////////////////////////////////////////////////////////////////////////
! generate the interpolation info for the im_traj_x with perturbation method
subroutine generate_perturb_x(ts)
    use mod_interp_x, only: n_samp, samp_dt, samp_t, samp_x, b, c, d
    implicit none
    double complex, intent(in):: ts
    double precision:: ti
    double complex:: c2
    double complex, external:: interp_perturb_vx, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i
    ti = dimag(ts)
    samp_dt = (ti - 0d0) / (n_samp - 1d0)
    c2 = - eye*dimag( simpson_sub( ti, 0d0, 1000, ts, interp_perturb_vx ) )
    forall(i = 1:n_samp) samp_t(i) = 0d0 + samp_dt * (i - 1d0)
    do i = 1, n_samp
        samp_x(i) = simpson_sub(ti, samp_t(i), 1000, ts, interp_perturb_vx) + c2
    end do
    call spline(n_samp, samp_t, samp_x, b, c, d)
    return
end subroutine generate_perturb_x

! ////////////////////////////////////////////////////////////////////////////////
! the final interpolation function for the im_traj_x with perturbation method
double complex function interp_perturb_x(t)
    use mod_interp_x, only: samp_dt, samp_t, samp_x, b, c, d
    implicit none
    double precision, intent(in):: t
    integer:: index
    index = ceiling( t / samp_dt )
    if(index .eq. 0) index = 1
    interp_perturb_x = samp_x(index) + b(index) * ( t - samp_t(index) ) + &
         c(index) * ( t - samp_t(index) )**2 + d(index) * ( t - samp_t(index) )**3
end function interp_perturb_x

! ////////////////////////////////////////////////////////////////////////////////
! plot the sub-barrier trajectories with the perturbation method
subroutine plot_sub_traj_ptb_1_interp(ts)
    implicit none
    double complex, intent(in):: ts
    integer:: i
    integer, parameter:: n_samp = 200
    double precision:: ti, dt, t
    double complex, external:: perturb_fx, interp_perturb_vx, interp_perturb_x

    ti = dimag(ts)-1e-6
    dt = ( ti - 0d0 ) / ( n_samp - 1d0 )

    open(101, file="dat/sub_traj_interp.dat")
    do i = 1, n_samp
        t = 0d0 + ( i - 1d0 ) * dt
        write(101, '(7(e15.8,1x))'), t, perturb_fx(t, ts),  &
             interp_perturb_vx(t, ts), interp_perturb_x(t)
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
    double complex:: c1, c2, vx(n_samp), x(n_samp)
    double complex, external:: perturb_fx, interp_perturb_vx, simpson_sub
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    integer:: i

    ti = dimag(ts)
    dt = (ti - 0d0) / (n_samp - 1d0)
    c1 = - eye * dimag( eye * simpson_sub( ti, 0d0, 1000, ts, perturb_fx ) )
    c2 = - eye * dimag( simpson_sub( ti, 0d0, 1000, ts, interp_perturb_vx ) )
    forall(i = 1:n_samp) t(i) = 0d0 + dt * (i - 1d0)
    do i = 1, n_samp
        vx(i) = eye * simpson_sub(ti, t(i), 1000, ts, perturb_fx) + c1
        x(i) = simpson_sub(ti, t(i), 1000, ts, interp_perturb_vx) + c2
    end do

    open(101, file="dat/sub_traj_samp.dat")
    do i = 1, n_samp
        write(101, '(5(e15.8,1x))'), t(i), vx(i), x(i)
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
