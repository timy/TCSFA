#include '../include/inc_rk4.h'
#include '../include/inc_atom.h'
#include '../include/inc_misc.h'
#include '../include/inc_field.h'

subroutine rk4_prop( t0, tp, x0, vx0, z0, vz0, ierr, w, px, pz, L, n_near_core, tag )

    implicit none
    double precision:: t0, tp, x0, vx0, z0, vz0, px, pz, L
    double complex:: w, w_near_core, dt
    double precision:: w_tail
    integer:: ierr
    integer:: tag

    integer, parameter:: nt = RK4_NT
    integer, parameter:: ne = RK4_NE
    double precision, parameter:: eps = RK4_EPS
    double precision, parameter:: charge = ATOM_CHARGE_Z
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: threshold = RK4_R_THRESHOLD

    double precision:: h, h0, h_old, t, t_old
    double precision:: y(ne), y0(ne), y_old(ne), y_mid(ne), y1(ne), y2(ne), y_new(ne)
    double precision:: ratio, r
    integer:: n_step, n_substep, n_near_core, b_near_core
    double complex, external:: action_w_re
    integer:: b_last
    
    
    
    h0 = ( tp - t0 ) / ( nt - 1d0 )
    y0(1) = x0; y0(2) = vx0; y0(3) = z0; y0(4) = vz0


#ifdef MISC_PLOT

    call rk4_plot_init(tag)
    call rk4_plot_open_file()
#endif
#if MISC_PRINT > 3
    print*, 'rk4_re(): Initialization of RK4'
    print*, 't0', t0, 'tp', tp;
    print*, 'x0', x0, 'vx0', vx0;
    print*, 'z0', z0, 'vz0', vz0;
    print*, 'h0', h0
#endif

    t = t0
    h = h0
    y = y0
    w = dcmplx( 0d0, 0d0 )
    n_step = 0
    n_substep = 0
    n_near_core = 0
    b_last = 0

10  n_step = n_step + 1
    t_old = t
    h_old = h
    y_old = y

    ! first try
    call rk4( h, t, y, y1 )
    
    ! second try
20  h = 0.5d0 * h
    if( h < eps ) then ! too small steps
        ierr = 2
        return
    end if
    t = t_old
    y = y_old
    call rk4( h, t, y, y_mid )
    t = t + h
    y = y_mid
    call rk4( h, t, y, y2 )

    ratio = maxval( dabs( y2 - y1 ) ) / eps
    
    if( ratio > 1 ) then

        y1 = y_mid
        n_substep = n_substep + 1
        go to 20

    else

        h = h * 2d0
        t = t_old + h
        y = y2
       
        w = w + action_W_re( h, t, y_old, y )


#ifdef MISC_PLOT
        call rk4_plot_write( n_step, t, y, h, n_substep )
#endif
        n_substep = 0

        ! adjust the step length to accelerate the calculation
        if( (ratio < 0.5d0) .and. (b_last .eq. 0) ) then
            h = h * 2d0
        end if

        ! if t > tp for the next step 
        if( t + h > tp ) then 
            h = tp - t
            b_last = 1
            ! if this step has been the end of the laser pulse
            if( dabs(h) < 1d-10 ) go to 30
        end if

        ! if enter the near-core region
        r = dsqrt( y(1)*y(1) + y(3)*y(3) )
        if( r < threshold ) then
!            call kepler_near_core( y(3), y(1), y(4), y(2), charge, Ip, dt, &
!                  y_new(3), y_new(1), y_new(4), y_new(2), w_near_core )
!            t = t + dreal( dt )
!            y = y_new
!            w = w + w_near_core
            if( b_near_core .eq. 0 ) then
                n_near_core = n_near_core + 1
                b_near_core = 1
            end if
        else
            b_near_core = 0
        end if


        
        go to 10

    end if

    

30  call coulomb_tail( y(3), y(1), y(4), y(2), charge, L, pz, px, w_tail, ierr )

    ! energy is less than 0: ierr == 3
    if( ierr > 0 ) then 
        return
    end if

    ! W_cen for the W-representation
    if( REPRESENTATION_OPT == 'W' ) then
        w = w + dcmplx( w_tail, 0d0 );
    end if

#ifdef MISC_PLOT
    call rk4_plot_close_file()
#endif
 
    ierr = 0

    return
    
end subroutine rk4_prop


! ----------------------------------------------------------------------
subroutine rk4( h, t0, y0, y )

    implicit none

    integer, parameter:: ne = RK4_NE
    double precision:: h, t0, y0(ne), y(ne)
    double precision:: t, k1(ne), k2(ne), k3(ne), k4(ne)

    t = t0
    y = y0

    call newton_equation_re( ne, t, y, k1 )

    t = t0 + 0.5d0 * h

    y = y0 + 0.5d0 * h * k1
    call newton_equation_re( ne, t, y, k2 )

    y = y0 + 0.5d0 * h * k2
    call newton_equation_re( ne, t, y, k3 )

    t = t0 + h
    y = y0 + h * k3
    call newton_equation_re( ne, t, y, k4 )

    y = y0 + h * ( k1 + 2d0 * ( k2 + k3 ) + k4 ) / 6d0

    return
end subroutine rk4

! ----------------------------------------------------------------------

double complex function action_W_re( h, t, y_old, y )
    implicit none
    double precision, parameter:: ne = RK4_NE
    double precision, parameter:: Ip = IONIZATION_IP
    double precision, parameter:: charge = ATOM_CHARGE_Z;
    double precision, intent(in):: h, t, y_old(ne), y(ne)
    double precision:: x, z, vx, vz, ax, az, Ex, Ez, v2, r, energy, t_mid
    double complex, external:: PULSE_E_X, PULSE_E_Z


    x = 0.5d0 * ( y_old(1) + y(1) )
    vx = 0.5d0 * ( y_old(2) + y(2) )
    z = 0.5d0 * ( y_old(3) + y(3) )
    vz = 0.5d0 * ( y_old(4) + y(4) )
    v2 = vz * vz + vx * vx
    r = dsqrt( x*x + z*z );

    
    if( REPRESENTATION_OPT == 'S' ) then
                        
        energy = 0.5d0 * v2 - charge / r;
        ! here modified for urgent use of screened potential !!!
        !energy = 0.5d0 * v2 - 1d0 / r - (charge - 1d0) * dexp(-5d0*r) / r

    elseif( REPRESENTATION_OPT == 'W' ) then
        
        t_mid = t - 0.5 * h;
        Ex = dreal( PULSE_E_X( dcmplx(t_mid) ) );
        Ez = dreal( PULSE_E_Z( dcmplx(t_mid) ) );
        ax = -Ex - charge * x / r**3;
        az = -Ez - charge * z / r**3;
        energy = 0.5d0 * v2 - charge / r + Ex * x + Ez * z + ax * x + az * z;

    end if
    
    action_W_re =  dcmplx( ( energy + Ip ) * h, 0d0 );

end function action_W_re
