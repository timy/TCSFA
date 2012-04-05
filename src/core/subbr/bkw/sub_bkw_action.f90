#include '../../../include/inc_atom.h'
#include '../../../include/inc_rk4.h'
#include '../../../include/inc_field.h'

double complex function action_w_sub_bkw( n_step, t, z, x, vz, vx ) result(w_sub)
    use mod_sub_bkw_prop, only: m_t0, m_ti
    implicit none
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    integer, intent(in):: n_step
    double precision, intent(in):: t(RK4_NT)
    double complex, intent(in):: z(RK4_NT), x(RK4_NT), vz(RK4_NT), vx(RK4_NT)
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double complex:: vzm, vxm, zm, xm
    integer:: i
    w_sub = dcmplx(0d0, 0d0)
    do i = 1, n_step
        vzm = 0.5d0*(vz(i)+vz(i+1))
        vxm = 0.5d0*(vx(i)+vx(i+1))
        zm = 0.5d0*(z(i)+z(i+1))
        xm = 0.5d0*(x(i)+x(i+1))
        w_sub = w_sub + (t(i+1) - t(i)) * (0.5*(vzm*vzm+vxm*vxm) - charge / cdsqrt(zm*zm+xm*xm) + Ip)
    end do
    w_sub = -w_sub ! due to the backward propagation
    w_sub = - eye * w_sub
end function action_w_sub_bkw
