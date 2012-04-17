#include '../../include/inc_field.h'

double complex function pulse_A_z_const( t ) result(A_z)
    implicit none;
    double complex, intent(in):: t;
    
    A_z = -( ( E0 * cdsin(OM*t) ) / ( OM * dsqrt( 1d0 + XI**2 ) ) )

    return;
end function pulse_A_z_const




double complex function pulse_A_x_const( t ) result(A_x)
    implicit none;
    double complex, intent(in):: t;

    A_x = - ( ( E0 * XI * cdcos( OM * t) ) / ( OM * dsqrt( 1 + XI**2 ) ) )
    
    return;    
end function pulse_A_x_const




double complex function pulse_E_z_const( t ) result(E_z)
    implicit none;
    double complex, intent(in):: t;

    E_z = ( E0 * cdcos(OM*t) ) / dsqrt( 1d0 + XI**2 )

    return;
end function pulse_E_z_const




double complex function pulse_E_x_const( t ) result(E_x)
    implicit none;
    double complex, intent(in):: t;
    
    E_x = - ( E0 * XI * cdsin(OM*t) ) / dsqrt( 1d0 + XI**2 )

    return;
end function pulse_E_x_const




double complex function pulse_alpha_z_const( t_ ) result(alpha_z)
    implicit none;
    double complex, intent(in):: t_;
    double complex:: t;

    if( dreal(t_) < 0 ) then

        alpha_z = (0d0, 0d0);
        return;
        
    else if( dreal(t_) > Tp ) then
        
        t = dcmplx(Tp, 0d0);
        
    else
        
        t = t_;
        
    end if
    
    alpha_z = (E0*(cdcos(OM*t) - dcos(OM*0d0)))/(OM**2*dsqrt(1 + XI**2));

    return;
end function pulse_alpha_z_const




double complex function pulse_alpha_x_const( t_ ) result(alpha_x)
    implicit none;
    double complex, intent(in):: t_;
    double complex:: t;
    
    if( dreal(t_) < 0 ) then
        
        alpha_x = (0d0, 0d0);
        return;
        
    else if( dreal(t_) > Tp ) then
        
        t = dcmplx(Tp, 0d0);
        
    else
    
        t = t_;
        
    end if
    alpha_x = -((E0*XI*(cdsin(OM*t) - dsin(OM*0d0)))/(OM**2*dsqrt(1 + XI**2)));
    
    return;
end function pulse_alpha_x_const


double complex function pulse_alpha2_z_const( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = -(E0**2*(-2*OM*t + cdsin(2*OM*t)))/(4.*OM**3)
    else if( t_ > Tp) then
        f = -(E0**2*(-2*OM*Tp + dsin(2*OM*Tp)))/(4.*OM**3)
    end if


    return;
end function pulse_alpha2_z_const




double complex function pulse_alpha2_x_const( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    f = dcmplx(0d0)
    return;
end function pulse_alpha2_x_const
