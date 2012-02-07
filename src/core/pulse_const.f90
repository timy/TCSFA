double complex function pulse_A_z_const( t ) result(A_z)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    
    A_z = -( ( E0 * cdsin(om*t) ) / ( om * dsqrt( 1d0 + xi**2 ) ) )

    return;
end function pulse_A_z_const




double complex function pulse_A_x_const( t ) result(A_x)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    A_x = - ( ( E0 * xi * cdcos( om * t) ) / ( om * dsqrt( 1 + xi**2 ) ) )
    
    return;    
end function pulse_A_x_const




double complex function pulse_E_z_const( t ) result(E_z)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    E_z = ( E0 * cdcos(om*t) ) / dsqrt( 1d0 + xi**2 )

    return;
end function pulse_E_z_const




double complex function pulse_E_x_const( t ) result(E_x)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    
    E_x = - ( E0 * xi * cdsin(om*t) ) / dsqrt( 1d0 + xi**2 )

    return;
end function pulse_E_x_const




double complex function pulse_alpha_z_const( t_ ) result(alpha_z)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t_;
    double complex:: t;

    if( dreal(t_) < 0 ) then

        alpha_z = (0d0, 0d0);
        return;
        
    else if( dreal(t_) > tp ) then
        
        t = dcmplx(tp, 0d0);
        
    else
        
        t = t_;
        
    end if
    
    alpha_z = (E0*(cdcos(om*t) - cdcos(om*t0)))/(om**2*dsqrt(1 + xi**2));

    return;
end function pulse_alpha_z_const




double complex function pulse_alpha_x_const( t_ ) result(alpha_x)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t_;
    double complex:: t;
    
    if( dreal(t_) < 0 ) then
        
        alpha_x = (0d0, 0d0);
        return;
        
    else if( dreal(t_) > tp ) then
        
        t = dcmplx(tp, 0d0);
        
    else
    
        t = t_;
        
    end if
    alpha_x = -((E0*xi*(cdsin(om*t) - cdsin(om*t0)))/(om**2*dsqrt(1 + xi**2)));
    
    return;
end function pulse_alpha_x_const


double complex function pulse_alpha2_z_const( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = -(E0**2*(-2*om*t + cdsin(2*om*t)))/(4.*om**3)
    else if( t_ > Tp) then
        f = -(E0**2*(-2*om*Tp + dsin(2*om*Tp)))/(4.*om**3)
    end if


    return;
end function pulse_alpha2_z_const




double complex function pulse_alpha2_x_const( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    f = dcmplx(0d0)
    return;
end function pulse_alpha2_x_const
