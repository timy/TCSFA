#include '../../../include/inc_atom.h'
#include '../../../include/inc_rk4.h'
#include '../../../include/inc_field.h'

! Note: x1_x2_x3 for the function name,
! x1: whether have 1/r term in the action
! x2: whether use 1st-order traj for 0.5 v^2
! x3: whether use 1st-order traj for 1/r

! W_sub with H = 0.5*(v0+v1)^2 - 1/|r0+r1|
double complex function action_w_sub_1_1_1( n_step, t, z1, x1, vz1, vx1, ts ) result(w_sub)
    use mod_sub_ptb_prop, only: m_t0, m_ti
    implicit none
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    integer, intent(in):: n_step
    double precision, intent(in):: t(RK4_NT)
    double complex, intent(in):: z1(RK4_NT), x1(RK4_NT), vz1(RK4_NT), vx1(RK4_NT)
    double complex, intent(in):: ts
    double complex, external:: sub_traj_z_0, sub_traj_x_0, sub_traj_vz_0, sub_traj_vx_0
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: t_mid, z, x, vz, vx
    double precision:: t0
    integer:: i
    t0 = dreal(ts)
    w_sub = dcmplx(0d0, 0d0)
    do i = 1, n_step
        t_mid = dcmplx( t0, 0.5d0 * ( t(i) + t(i+1) ) )
        vz = sub_traj_vz_0(t_mid) + 0.5d0*(vz1(i)+vz1(i+1))
        vx = sub_traj_vx_0(t_mid) + 0.5d0*(vx1(i)+vx1(i+1))
        z = sub_traj_z_0(t_mid, ts) + 0.5d0*(z1(i)+z1(i+1))
        x = sub_traj_x_0(t_mid, ts) + 0.5d0*(x1(i)+x1(i+1))
        w_sub = w_sub + (t(i+1) - t(i)) * (0.5*(vz*vz+vx*vx) - charge / cdsqrt(z*z+x*x) + Ip)
    end do
    w_sub = - eye * w_sub
end function action_w_sub_1_1_1

! W_sub with H = 0.5*v0^2; there is also an analytical version
double complex function action_w_sub_0_0_0( n_step, t, z1, x1, vz1, vx1, ts ) result(w_sub)
    use mod_sub_ptb_prop, only: m_t0, m_ti
    implicit none
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    integer, intent(in):: n_step
    double precision, intent(in):: t(RK4_NT)
    double complex, intent(in):: z1(RK4_NT), x1(RK4_NT), vz1(RK4_NT), vx1(RK4_NT)
    double complex, intent(in):: ts
    double complex, external:: sub_traj_z_0, sub_traj_x_0, sub_traj_vz_0, sub_traj_vx_0
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: t_mid, z, x, vz, vx
    double precision:: t0
    integer:: i
    t0 = dreal(ts)
    w_sub = dcmplx(0d0, 0d0)
    do i = 1, n_step
        t_mid = dcmplx( t0, 0.5d0 * ( t(i) + t(i+1) ) )
        vz = sub_traj_vz_0(t_mid)
        vx = sub_traj_vx_0(t_mid)
        z = sub_traj_z_0(t_mid, ts)
        x = sub_traj_x_0(t_mid, ts)
        w_sub = w_sub + (t(i+1) - t(i)) * (0.5*(vz*vz+vx*vx) + Ip)
    end do
    w_sub = - eye * w_sub
end function action_w_sub_0_0_0

! W_sub with H = 0.5*v0^2 - 1/|r0|; there is also an analytical version
double complex function action_w_sub_1_0_0( n_step, t, z1, x1, vz1, vx1, ts ) result(w_sub)
    use mod_sub_ptb_prop, only: m_t0, m_ti
    implicit none
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    integer, intent(in):: n_step
    double precision, intent(in):: t(RK4_NT)
    double complex, intent(in):: z1(RK4_NT), x1(RK4_NT), vz1(RK4_NT), vx1(RK4_NT)
    double complex, intent(in):: ts
    double complex, external:: sub_traj_z_0, sub_traj_x_0, sub_traj_vz_0, sub_traj_vx_0
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: t_mid, z, x, vz, vx
    double precision:: t0
    integer:: i
    t0 = dreal(ts)
    w_sub = dcmplx(0d0, 0d0)
    do i = 1, n_step
        t_mid = dcmplx( t0, 0.5d0 * ( t(i) + t(i+1) ) )
        vz = sub_traj_vz_0(t_mid)
        vx = sub_traj_vx_0(t_mid)
        z = sub_traj_z_0(t_mid, ts)
        x = sub_traj_x_0(t_mid, ts)
        w_sub = w_sub + (t(i+1) - t(i)) * (0.5*(vz*vz+vx*vx) - charge / cdsqrt(z*z+x*x) + Ip)
    end do
    w_sub = - eye * w_sub
end function action_w_sub_1_0_0

! W_sub with H = 0.5*v0^2 - 1/|r0+r1|
double complex function action_w_sub_1_0_1( n_step, t, z1, x1, vz1, vx1, ts ) result(w_sub)
    use mod_sub_ptb_prop, only: m_t0, m_ti
    implicit none
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    integer, intent(in):: n_step
    double precision, intent(in):: t(RK4_NT)
    double complex, intent(in):: z1(RK4_NT), x1(RK4_NT), vz1(RK4_NT), vx1(RK4_NT)
    double complex, intent(in):: ts
    double complex, external:: sub_traj_z_0, sub_traj_x_0, sub_traj_vz_0, sub_traj_vx_0
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: t_mid, z, x, vz, vx
    double precision:: t0
    integer:: i
    t0 = dreal(ts)
    w_sub = dcmplx(0d0, 0d0)
    do i = 1, n_step
        t_mid = dcmplx( t0, 0.5d0 * ( t(i) + t(i+1) ) )
        vz = sub_traj_vz_0(t_mid)
        vx = sub_traj_vx_0(t_mid)
        z = sub_traj_z_0(t_mid, ts) + 0.5d0*(z1(i)+z1(i+1))
        x = sub_traj_x_0(t_mid, ts) + 0.5d0*(x1(i)+x1(i+1))
        w_sub = w_sub + (t(i+1) - t(i)) * (0.5*(vz*vz+vx*vx) - charge / cdsqrt(z*z+x*x) + Ip)
    end do
    w_sub = - eye * w_sub
end function action_w_sub_1_0_1
