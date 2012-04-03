#include '../../include/inc_rk4.h'
#include '../../include/inc_atom.h'
#include '../../include/inc_misc.h'
#include '../../include/inc_field.h'

module mod_sub_prop
    implicit none
    double precision:: m_t0, m_ti
end module mod_sub_prop

subroutine sub_prop( ts, ierr, w, tag )
    use mod_sub_prop, only: m_t0, m_ti
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
    double complex:: y_actual(ne), y_old_actual(ne)
    double precision:: ratio
    integer:: n_step, n_substep
    double complex, external:: action_w_sub
    integer:: b_last
    double complex, external:: sub_traj_x_0, sub_traj_z_0
    double complex, external:: sub_traj_vx_0, sub_traj_vz_0
    
    m_t0 = dreal(ts)
    m_ti = dimag(ts)

    h0 = (0d0 - m_ti) / (nt - 1d0)
    y0(1) = sub_traj_x_0( ts, ts )
    y0(2) = eye * sub_traj_vx_0( ts )
    y0(3) = sub_traj_z_0( ts, ts )
    y0(4) = eye * sub_traj_vz_0( ts )

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
        y_old_actual(1) = y_old(1); y_old_actual(3) = y_old(3)
        y_old_actual(2) = -eye*y_old(2); y_old_actual(4) = -eye*y_old(4)
        y_actual(1) = y(1); y_actual(3) = y(3)
        y_actual(2) = -eye*y(2); y_actual(4) = -eye*y(4)
        w = w + action_w_sub( h, t, y_old_actual, y_actual )
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
    return
    
end subroutine sub_prop

! ----------------------------------------------------------------------
double complex function action_w_sub( h, t, y_old, y )
    use mod_sub_prop, only: m_t0, m_ti
    implicit none
    double precision, parameter:: ne = RK4_NE
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    double precision, intent(in):: h, t
    double complex, intent(in):: y_old(ne), y(ne)
    double precision:: t_mid
    double complex:: x, z, vx, vz, ax, az, Ex, Ez, v2, r, energy
    double complex, external:: PULSE_E_X, PULSE_E_Z


    x = 0.5d0 * ( y_old(1) + y(1) )
    vx = 0.5d0 * ( y_old(2) + y(2) )
    z = 0.5d0 * ( y_old(3) + y(3) )
    vz = 0.5d0 * ( y_old(4) + y(4) )
    v2 = vz * vz + vx * vx
    r = cdsqrt( x*x + z*z );

    
    if( REPRESENTATION_OPT == 'S' ) then
                        
        energy = 0.5d0 * v2 - charge / r;
        ! here modified for urgent use of screened potential !!!
        !energy = 0.5d0 * v2 - 1d0 / r - (charge - 1d0) * dexp(-5d0*r) / r

    elseif( REPRESENTATION_OPT == 'W' ) then
        
        t_mid = t - 0.5 * h;
        Ex = PULSE_E_X( dcmplx(m_t0, t_mid) );
        Ez = PULSE_E_Z( dcmplx(m_t0, t_mid) );
        ax = -Ex - charge * x / r**3;
        az = -Ez - charge * z / r**3;
        energy = 0.5d0 * v2 - charge / r + Ex * x + Ez * z + ax * x + az * z;

    end if
    
    action_W_sub =  - ( energy + Ip ) * h; ! whether there is some "eye" here?

end function action_w_sub

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
    use mod_sub_prop, only: m_t0
    implicit none
    integer, intent(in):: ne
    double precision, intent(in):: t
    double complex, intent(in):: y(ne)
    double complex, intent(out):: dy(ne)
    double complex, external:: PULSE_E_Z, PULSE_E_X
    double complex:: r3
    double precision:: charge

    charge = ATOM_CHARGE_Z
    r3 = ( y(1)*y(1) + y(3)*y(3) ) ** 1.5
    ! y(1) = x; y(2) = vx; y(3) = z; y(4) = vz;
    ! EOM. Plus signs, for the imaginary integration path.
    dy(1) = y(2)
    dy(2) = PULSE_E_X( dcmplx(m_t0, t) ) + y(1) * charge / r3
    dy(3) = y(4)
    dy(4) = PULSE_E_Z( dcmplx(m_t0, t) ) + y(3) * charge / r3
    return
end subroutine newton_equation_sub
