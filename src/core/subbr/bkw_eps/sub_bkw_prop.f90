#include '../../../include/inc_rk4.h'
#include '../../../include/inc_atom.h'
#include '../../../include/inc_misc.h'
#include '../../../include/inc_field.h'

module mod_sub_bkw_prop_eps
    implicit none
    double precision:: m_t0, m_ti
end module mod_sub_bkw_prop_eps

subroutine sub_bkw_prop_eps( ts, ierr, z_t0, x_t0, vz_t0, vx_t0, w, err_spe, tag )
    use mod_sub_bkw_prop_eps, only: m_t0, m_ti
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
    double complex, external:: action_w_sub_bkw_eps
    integer:: b_last
    double complex, external:: sub_traj_x_0, sub_traj_z_0
    double complex, external:: sub_traj_vx_0, sub_traj_vz_0
    double precision:: array_t(nt)
    double complex:: array_z(nt), array_x(nt), array_vz(nt), array_vx(nt)
    double complex, intent(out):: z_t0, x_t0, vz_t0, vx_t0
    double precision, intent(out):: err_spe
    external:: newton_equation_sub_bkw_eps

    integer, parameter:: n_interval = 10
    integer:: i_interval
    double precision:: h_interval, interval_ta, interval_tb

    m_t0 = dreal(ts)
    m_ti = dimag(ts)

    h_interval = ( m_ti - 0d0 ) / n_interval
    w = dcmplx( 0d0, 0d0 )

#ifdef MISC_PLOT_TRAJ
    call rk4_plot_init(1, tag)
    call rk4_plot_open_file()
#endif

    do i_interval = 1, n_interval
        
        interval_ta = h_interval * ( i_interval - 1d0 )
        interval_tb = interval_ta + h_interval 
        h0 = (interval_tb - interval_ta) / (nt - 1d0)

        array_t(1) = interval_ta
        array_x(1) = sub_traj_x_0(dcmplx(m_t0, interval_ta), ts)
        array_z(1) = sub_traj_z_0(dcmplx(m_t0, interval_ta), ts)
        array_vx(1) = sub_traj_vx_0(dcmplx(m_t0, interval_ta))
        array_vz(1) = sub_traj_vz_0(dcmplx(m_t0, interval_ta))

        y0(1) = array_x(1)
        y0(2) = eye * array_vx(1)
        y0(3) = array_z(1)
        y0(4) = eye * array_vz(1)

#if MISC_PRINT > 10
        write(*,'(a,x,a)'), 'Print in『sub_prop』', repeat('-', 40)
        write(*, '(a,2(e15.8,x))'), 'h0:', h0
        write(*, '(a,2(e15.8,x))'), 'z~(ti):', y0(3)
        write(*, '(a,2(e15.8,x))'), 'x~(ti):', y0(1)
        write(*, '(a,2(e15.8,x))'), 'vz~(ti):', y0(4)
        write(*, '(a,2(e15.8,x))'), 'vx~(ti):', y0(2)
#endif

        t = interval_ta
        h = h0
        y = y0
        n_step = 0
        n_substep = 0
        

10      n_step = n_step + 1
        t_old = t
        h_old = h
        y_old = y
    
        ! first try
        call rk4_sub( h, t, y, y1, newton_equation_sub_bkw_eps )
    
    ! second try
20      b_last = 0
        h = 0.5d0 * h
        if( dabs(h) < eps ) then ! too small steps
            ierr = 2
            return
        end if
        t = t_old
        y = y_old
        call rk4_sub( h, t, y, y_mid, newton_equation_sub_bkw_eps )
        t = t + h
        y = y_mid
        call rk4_sub( h, t, y, y2, newton_equation_sub_bkw_eps )

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
            array_x(n_step+1) = y_actual(1)
            array_z(n_step+1) = y_actual(3)
            array_vx(n_step+1) = y_actual(2)
            array_vz(n_step+1) = y_actual(4)

#ifdef MISC_PLOT_TRAJ
            call rk4_plot_write_cmplx( n_step, t, y_actual, h, n_substep )
#endif
            n_substep = 0

            ! adjust the step length to accelerate the calculation
            if( (ratio < 0.5d0) .and. (b_last .eq. 0) ) then
                h = h * 2d0
            end if

            ! if t > tp for the next step 
            if( t + h > interval_tb ) then 
                h = interval_tb - t
                b_last = 1
                ! if this step has been the end of the laser pulse
                if( dabs(h) < 1d-10 ) go to 30
            end if

            go to 10 ! goto the next step

        end if

30      ierr = 0
        w = w + action_w_sub_bkw_eps( n_step, array_t, array_z, array_x, array_vz, array_vx )
    end do

#ifdef MISC_PLOT_TRAJ
    call rk4_plot_close_file()
#endif

    z_t0 = sub_traj_z_0( dcmplx(m_t0, 0d0), ts )
    x_t0 = sub_traj_x_0( dcmplx(m_t0, 0d0), ts )
    vz_t0 = sub_traj_vz_0( dcmplx(m_t0, 0d0) )
    vx_t0 = sub_traj_vx_0( dcmplx(m_t0, 0d0) )
    err_spe = 0d0
    return    
end subroutine sub_bkw_prop_eps

! --------------------------------------------------------------------------------
subroutine newton_equation_sub_bkw_eps( ne, t, y, dy )
    use mod_sub_bkw_prop_eps, only: m_t0,  m_ti
    implicit none
    integer, intent(in):: ne
    double precision, intent(in):: t
    double complex, intent(in):: y(ne)
    double complex, intent(out):: dy(ne)
    double complex:: r3
    double precision:: charge
    double complex, external:: PULSE_E_Z;

    charge = ATOM_CHARGE_Z
    r3 = ( y(1)*y(1) + y(3)*y(3) ) ** 1.5
    ! y(1) = x; y(2) = vx; y(3) = z; y(4) = vz;
    ! EOM. Plus signs, for the imaginary integration path.
    dy(1) = y(2)
    dy(2) = charge * y(1) / r3
    dy(3) = y(4)
    dy(4) = PULSE_E_Z(dcmplx(m_t0, t)) + charge * y(3) / r3
    return
end subroutine newton_equation_sub_bkw_eps
