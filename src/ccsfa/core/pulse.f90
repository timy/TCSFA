#define PULSE_NT 200
#define PULSE_FILE_ID 101
#define PULSE_FILE_NAME 'dat/pulse.dat'


module mod_pulse

    implicit none;

    double precision:: E0;
    double precision:: om;
    double precision:: nc;
    double precision:: xi;
    double precision, parameter::  PI = 2d0*asin(1d0);

    double precision:: tp;

    ! t0 is used only when calculating alpha, actually it has no influence whatever value
    ! it has, since it vanishes in the calculation of x0 and z0
    double complex:: t0 = ( 0d0, 0d0 );     

contains

    subroutine set_pulse( E0_, om_, nc_, xi_ )

        implicit none;
        double precision, intent(in):: E0_, om_, nc_, xi_

        E0 = E0_;
        om = om_;
        nc = nc_;
        xi = xi_;
        tp = 2d0*PI*nc / om;

        return;
    end subroutine set_pulse


    subroutine set_pulse_t0( t0_ )

        implicit none;
        double complex, intent(in):: t0_

        t0 = t0_;

        return;
    end subroutine set_pulse_t0


end module mod_pulse




double complex function pulse_A_z( t )
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        pulse_A_z = ( 0d0, 0d0 )
    else
        pulse_A_z = -((E0*cdsin(om*t))/(om*dsqrt(1 + xi**2)));
    end if

    return;
end function pulse_A_z




double complex function pulse_A_x( t )
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        pulse_A_x = ( 0d0, 0d0 )
    else
        pulse_A_x = -((E0*xi*cdcos(om*t))/(om*dsqrt(1 + xi**2)));
    end if

    return;    
end function pulse_A_x




double complex function pulse_E_z( t )
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;


    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        pulse_E_z = ( 0d0, 0d0 )
    else
        pulse_E_z = (E0*cdcos(om*t))/dsqrt(1 + xi**2);
    end if

    return;
end function pulse_E_z




double complex function pulse_E_x( t )
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;


    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        pulse_E_x = ( 0d0, 0d0 )
    else
        pulse_E_x = -((E0*xi*cdsin(om*t))/dsqrt(1 + xi**2));
    end if

    return;
end function pulse_E_x




double complex function alpha_z( t_ )
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t_;
    double complex:: t;

    if( dreal(t_) < 0 ) then

        alpha_z = ( 0d0, 0d0 );
        return;

    else if( dreal(t_) > tp ) then

        t = dcmplx(tp, 0d0);

    else

        t = t_;

    end if

    alpha_z = (E0*(cdcos(om*t) - cdcos(om*t0)))/(om**2*dsqrt(1 + xi**2));

    return;
end function alpha_z




double complex function alpha_x( t_ )
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t_;
    double complex:: t;

    if( dreal(t_) < 0 ) then

        alpha_x = ( 0d0, 0d0 );
        return;

    else if( dreal(t_) > tp ) then

        t = dcmplx(tp, 0d0);

    else

        t = t_;

    end if
    
    alpha_x = -((E0*xi*(cdsin(om*t) - cdsin(om*t0)))/(om**2*dsqrt(1 + xi**2)));

    return;
end function alpha_x



double precision function get_tp()
    use mod_pulse, only: tp;
    implicit none;

    get_tp = tp;

    return;
end function get_tp




subroutine plot_pulse()
    use mod_pulse

    implicit none
    double complex, external:: pulse_A_z, pulse_A_x, pulse_E_z, pulse_E_x
    double precision:: dt
    integer:: i
    double complex:: t

    open( PULSE_FILE_ID, file=PULSE_FILE_NAME );
    dt = tp / ( PULSE_NT - 1d0 );
    
    do i = 1, PULSE_NT

        t = dcmplx( (i - 1d0) * dt, 0d0 );
        write( PULSE_FILE_ID, '(5(e15.8,2x))'), dreal(t), &
              dreal(pulse_A_z(t)), dreal(pulse_A_x(t)), &
              dreal(pulse_E_z(t)), dreal(pulse_E_x(t));
        
    end do

    close( PULSE_FILE_ID );

end subroutine plot_pulse
