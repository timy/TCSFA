#include '../../include/inc_rk4.h'
#include '../../include/inc_atom.h'
#include '../../include/inc_misc.h'
#include '../../include/inc_field.h'

#define ACTION_W_SUB action_w_sub_1_0_1


module mod_sub_ptb_prop
    implicit none
    double precision:: m_t0, m_ti
end module mod_sub_ptb_prop

subroutine sub_ptb_prop( ts, ierr, z_t0, x_t0, vz_t0, vx_t0, w, err_spe, tag )
    use mod_sub_ptb_prop, only: m_t0, m_ti
    implicit none
    double complex, intent(in):: ts
    double complex:: w
    integer:: ierr, tag

    integer, parameter:: nt = RK4_NT
    integer, parameter:: ne = RK4_NE
    double precision, parameter:: eps = RK4_EPS
    double precision, parameter:: charge = ATOM_CHARGE_Z
    double precision, parameter:: Ip = IONIZATION_IP
    double complex, parameter:: eye = dcmplx(0d0, 1d0)
    double precision:: h, h0, h_old, t, t_old
    double complex:: y(ne), y0(ne), y_old(ne), y_mid(ne), y1(ne), y2(ne)
    double complex:: y_actual(ne)
    double precision:: ratio
    integer:: n_step, n_substep
    double complex, external:: ACTION_W_SUB
    integer:: b_last
    double complex, external:: sub_traj_x_0, sub_traj_z_0
    double complex, external:: sub_traj_vx_0, sub_traj_vz_0
    double complex, external:: sub_ptb1_fz, sub_ptb1_fx
    double complex, external:: simpson_sub
    double precision:: array_t(nt)
    double complex:: array_z(nt), array_x(nt), array_vz(nt), array_vx(nt)
    double complex:: Bpz, Bpx, Cpz, Cpx
    double complex, intent(out):: z_t0, x_t0, vz_t0, vx_t0
    double precision, intent(out):: err_spe
    double complex:: vz_ts, vx_ts
    integer:: i

    m_t0 = dreal(ts)
    m_ti = dimag(ts)

    Bpz = -eye * simpson_sub(0d0, m_ti, 1000, ts, sub_ptb1_fz)
    Bpx = -eye * simpson_sub(0d0, m_ti, 1000, ts, sub_ptb1_fx)

    h0 = (0d0 - m_ti) / (nt - 1d0)
    y0(1) = dcmplx(0d0, 0d0)
    y0(2) = dcmplx(0d0, 0d0)
    y0(3) = dcmplx(0d0, 0d0)
    y0(4) = dcmplx(0d0, 0d0)

    array_t(1) = m_ti
    array_x(1) = y0(1)
    array_z(1) = y0(3)
    array_vx(1) = -eye*y0(2) - Bpx
    array_vz(1) = -eye*y0(4) - Bpz

#if MISC_PRINT > 0
    write(*,'(a,x,a)'), 'Print in『sub_prop』', repeat('-', 40)
    write(*, '(a,2(e15.8,x))'), 'h0:', h0
    write(*, '(a,2(e15.8,x))'), 'z~(ti):', y0(3)
    write(*, '(a,2(e15.8,x))'), 'x~(ti):', y0(1)
    write(*, '(a,2(e15.8,x))'), 'vz~(ti):', y0(4)
    write(*, '(a,2(e15.8,x))'), 'vx~(ti):', y0(2)
#endif

#ifdef MISC_PLOT_TRAJ
    call rk4_plot_init(1, tag)
    call rk4_plot_open_file()
#endif

    t = m_ti
    h = h0
    y = y0
    w = dcmplx( 0d0, 0d0 )
    n_step = 0
    n_substep = 0
    b_last = 0

10  n_step = n_step + 1
    t_old = t
    h_old = h
    y_old = y
    
    ! first try
    call rk4_sub( h, t, y, y1 )
    
    ! second try
20  h = 0.5d0 * h
    if( dabs(h) < eps ) then ! too small steps
        ierr = 2
        return
    end if
    t = t_old
    y = y_old
    call rk4_sub( h, t, y, y_mid )
    t = t + h
    y = y_mid
    call rk4_sub( h, t, y, y2 )

    ratio = maxval( cdabs( y2 - y1 ) ) / eps
    
    if( ratio > 1 ) then ! the error is larger than tolerance
        y1 = y_mid
        n_substep = n_substep + 1
        go to 20 ! halve the step length and continue
    else ! accurate enough
        h = h * 2d0
        t = t_old + h
        y = y2
       
        ! pay attention to the complex variables
 !       y_old_actual(1) = y_old(1); y_old_actual(3) = y_old(3)
 !       y_old_actual(2) = -eye*y_old(2); y_old_actual(4) = -eye*y_old(4)
        y_actual(1) = y(1); y_actual(3) = y(3)
        y_actual(2) = -eye*y(2); y_actual(4) = -eye*y(4)

        if(n_step > nt) then
            ierr = 10  ! too many steps to store for the sub-prop
            return
        end if

        array_t(n_step+1) = t
        array_x(n_step+1) = y(1)
        array_z(n_step+1) = y(3)
        array_vx(n_step+1) = -eye*y(2) - Bpx
        array_vz(n_step+1) = -eye*y(4) - Bpz

#ifdef MISC_PLOT_TRAJ
        call rk4_plot_write_cmplx( n_step, t, y_actual, h, n_substep )
#endif
        n_substep = 0

        ! adjust the step length to accelerate the calculation
        if( (ratio < 0.5d0) .and. (b_last .eq. 0) ) then
            h = h * 2d0
        end if

        ! if t > tp for the next step 
        if( t + h < 0d0 ) then 
            h = - t
            b_last = 1
            ! if this step has been the end of the laser pulse
            if( dabs(h) < 1d-10 ) go to 30
        end if

        go to 10 ! goto the next step

    end if

30    ierr = 0
#ifdef MISC_PLOT_TRAJ
    call rk4_plot_close_file()
#endif

    Cpz = dcmplx(0d0, 0d0)
    Cpx = dcmplx(0d0, 0d0)
    do i = 1, n_step
        Cpz = Cpz + ( array_vz(i+1) + array_vz(i) ) * ( array_t(i+1) - array_t(i) )
        Cpx = Cpx + ( array_vx(i+1) + array_vx(i) ) * ( array_t(i+1) - array_t(i) )
    end do
    Cpz = eye * ( - dimag( eye * 0.5d0 * Cpz ) + Bpz * m_ti )
    Cpx = eye * ( - dimag( eye * 0.5d0 * Cpx ) + Bpx * m_ti )

    forall(i=1:n_step+1) array_z(i) = array_z(i) - eye * Bpz * array_t(i) + Cpz
    forall(i=1:n_step+1) array_x(i) = array_x(i) - eye * Bpx * array_t(i) + Cpx

#ifdef MISC_PLOT_TRAJ
    call rk4_plot_init(1, tag)
    call rk4_plot_open_file()
    !
    do i = 1, n_step+1
        y_actual(1) = sub_traj_x_0(dcmplx(m_t0, array_t(i)), ts) + array_x(i)
        y_actual(2) = sub_traj_vx_0(dcmplx(m_t0, array_t(i))) + array_vx(i)
        y_actual(3) = sub_traj_z_0(dcmplx(m_t0, array_t(i)), ts) + array_z(i)
        y_actual(4) = sub_traj_vz_0(dcmplx(m_t0, array_t(i))) + array_vz(i)
        call rk4_plot_write_cmplx( i, array_t(i), y_actual, 0d0, 0 )
    end do
    call rk4_plot_close_file()
#endif

    z_t0 = sub_traj_z_0(dcmplx(m_t0, 0d0), ts)! + array_z(n_step+1)
    x_t0 = sub_traj_x_0(dcmplx(m_t0, 0d0), ts)! + array_x(n_step+1)
    vz_t0 = sub_traj_vz_0(dcmplx(m_t0, 0d0)) !+ array_vz(n_step+1)
    vx_t0 = sub_traj_vx_0(dcmplx(m_t0, 0d0)) !+ array_vx(n_step+1)
    w =  ACTION_W_SUB( n_step, array_t, array_z, array_x, array_vz, array_vx, ts )

    vz_ts = sub_traj_vz_0(ts) + array_vz(1)
    vx_ts = sub_traj_vx_0(ts) + array_vx(1)
    err_spe = cdabs(0.5 * (vz_ts*vz_ts + vx_ts*vx_ts) + Ip)
    return
    
end subroutine sub_ptb_prop

! --------------------------------------------------------------------------------
subroutine rk4_sub( h, t0, y0, y )

    implicit none
    integer, parameter:: ne = RK4_NE
    double precision:: h, t0, t
    double complex:: y0(ne), y(ne)
    double complex:: k1(ne), k2(ne), k3(ne), k4(ne)

    t = t0
    y = y0

    call newton_equation_sub( ne, t, y, k1 )

    t = t0 + 0.5d0 * h

    y = y0 + 0.5d0 * h * k1
    call newton_equation_sub( ne, t, y, k2 )

    y = y0 + 0.5d0 * h * k2
    call newton_equation_sub( ne, t, y, k3 )

    t = t0 + h
    y = y0 + h * k3
    call newton_equation_sub( ne, t, y, k4 )

    y = y0 + h * ( k1 + 2d0 * ( k2 + k3 ) + k4 ) / 6d0

    return
end subroutine rk4_sub

! --------------------------------------------------------------------------------
subroutine newton_equation_sub( ne, t, y, dy )
    use mod_sub_ptb_prop, only: m_t0,  m_ti
    implicit none
    integer, intent(in):: ne
    double precision, intent(in):: t
    double complex, intent(in):: y(ne)
    double complex, intent(out):: dy(ne)
    double complex:: r3
    double precision:: charge
    double complex:: ts
    double complex, external:: sub_ptb1_fz, sub_ptb1_fx

    ts = dcmplx( m_t0, m_ti )
    charge = ATOM_CHARGE_Z
    r3 = ( y(1)*y(1) + y(3)*y(3) ) ** 1.5
    ! y(1) = x; y(2) = vx; y(3) = z; y(4) = vz;
    ! EOM. Plus signs, for the imaginary integration path.
    dy(1) = y(2)
    dy(2) = - sub_ptb1_fx(t, ts)
    dy(3) = y(4)
    dy(4) = - sub_ptb1_fz(t, ts)
    return
end subroutine newton_equation_sub

! --------------------------------------------------------------------------------
double complex function sub_ptb1_fz( t, ts )
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex:: tc, r3, z0, x0
    double complex, external:: sub_traj_z_0, sub_traj_x_0
    
    tc = dcmplx( dreal(ts), t )
    z0 = sub_traj_z_0( tc, ts )
    x0 = sub_traj_x_0( tc, ts )
    
    r3 = (z0*z0 + x0*x0)**1.5
    sub_ptb1_fz = - ATOM_CHARGE_Z * z0 / r3 
    return
end function sub_ptb1_fz

double complex function sub_ptb1_fx( t, ts )
    implicit none
    double precision, intent(in):: t
    double complex, intent(in):: ts
    double complex:: tc, r3, z0, x0
    double complex, external:: sub_traj_z_0, sub_traj_x_0
    
    tc = dcmplx( dreal(ts), t )
    z0 = sub_traj_z_0( tc, ts )
    x0 = sub_traj_x_0( tc, ts )
    
    r3 = (z0*z0 + x0*x0)**1.5
    sub_ptb1_fx = - ATOM_CHARGE_Z * x0 / r3 
    return
end function sub_ptb1_fx
