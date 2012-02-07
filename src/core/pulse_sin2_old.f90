double complex function pulse_A_z_sin2( t ) result(A_z)
    use mod_pulse;
    implicit none;
    double complex, intent(in):: t;

    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        A_z = ( 0d0, 0d0 )
    else
!        pulse_A_z = -((E0*cdsin(om*t)*cdsin((om*t)/(2.*nc))**2)/(om*dsqrt(1 + xi**2)));
        A_z = -((E0*cdsin((om*t)/(2d0*nc))**2*cdsin(ph + om*t))/(om*dsqrt(1 + xi**2)));
    end if

    return;
end function pulse_A_z_sin2




double complex function pulse_A_x_sin2( t ) result(A_x)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        A_x = ( 0d0, 0d0 )
    else
        ! without phi
!        pulse_A_x = -((E0*xi*cdcos(om*t)*cdsin((-Pi/2. + om*t)/(2.*nc))**2)/ &
!              (om*dsqrt(1 + xi**2)))
        ! with phi
        A_x = -((E0*xi*cdcos(ph + om*t)*cdsin((-Pi/2. + om*t)/(2.*nc))**2)/ &
              (om*dsqrt(1 + xi**2)))
    end if
    
    return;    
end function pulse_A_x_sin2




double complex function pulse_E_z_sin2( t ) result(E_z)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        E_z = ( 0d0, 0d0 )
    else
        ! without phi
!        pulse_E_z = (E0*cdsin((om*t)/(2.*nc))*(cdcos((om*t)/(2.*nc))*cdsin(om*t) + &
!              nc*cdcos(om*t)*cdsin((om*t)/(2.*nc))))/(nc*dsqrt(1 + xi**2));
        ! with phi
        E_z = (E0*cdcos(ph + om*t)*cdsin((om*t)/(2.*nc))**2)/dsqrt(1 + xi**2) + &
              (E0*cdcos((om*t)/(2.*nc))*cdsin((om*t)/(2.*nc))*cdsin(ph + om*t))/(nc*dsqrt(1 + xi**2));
    end if

    return;
end function pulse_E_z_sin2




double complex function pulse_E_x_sin2( t ) result(E_x)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;

    if( dreal(t) < 0d0 .or. dreal(t) > tp ) then
        E_x = ( 0d0, 0d0 )
    else
        ! without phi
!        pulse_E_x = -((E0*xi*cdsin((Pi - 2*om*t)/(4.*nc))* &
!              (cdcos(om*t)*cdcos((Pi - 2*om*t)/(4.*nc)) + &
!              nc*cdsin(om*t)*cdsin((Pi - 2*om*t)/(4.*nc))))/(nc*dSqrt(1 + xi**2)));
        ! with phi
        E_x = (E0*xi*cdcos(ph + om*t)*cdcos((-Pi/2. + om*t)/(2.*nc))* &
              cdsin((-Pi/2. + om*t)/(2.*nc)))/(nc*dsqrt(1 + xi**2)) - &
              (E0*xi*cdsin(ph + om*t)*cdsin((-Pi/2. + om*t)/(2.*nc))**2)/dsqrt(1 + xi**2);
    end if

    return;
end function pulse_E_x_sin2




double complex function pulse_alpha_z_sin2( t_ ) result(alpha_z)
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


    if( dabs(nc - 1d0) < 1d-6 ) then

!        alpha_z = (E0*(4*cdcos(om*t) - cdcos(2*om*t) - 4*cdcos(om*t0) + &
!              cdcos(2*om*t0)))/(8.*om**2*dsqrt(1 + xi**2));

        alpha_z = (E0*(4*cdcos(ph + om*t) - cdcos(ph + 2*om*t) - 4*cdcos(ph + om*t0) &
              + cdcos(ph + 2*om*t0) + 2*om*t*dsin(ph) - 2*om*t0*dsin(ph)))/(8.*om**2*dsqrt(1 + xi**2))
    else

!        alpha_z = -(E0*(cdcos(om*t)*(1 - nc**2 + nc**2*cdcos((om*t)/nc)) + &
!              cdcos(om*t0)*(-1 + nc**2 - nc**2*cdcos((om*t0)/nc)) + &
!              nc*(cdsin(om*t)*cdsin((om*t)/nc) - cdsin(om*t0)*cdsin((om*t0)/nc))))/&
!              (2.*(-1 + nc**2)*om**2*dsqrt(1 + xi**2))

        alpha_z = -(E0*(-2*(-1 + nc**2)*cdcos(ph + om*t) + nc*(1 + nc) * cdcos(ph + &
              ((-1 + nc)*om*t)/nc) + (-1 + nc)*nc*cdcos(ph + ((1 + nc)*om*t)/nc) + &
              2*(-1 + nc**2)*cdcos(ph + om*t0) - nc*(1 + nc)*cdcos(ph + ((-1 + nc)*om*t0)/nc) &
              - (-1 + nc)*nc*cdcos(ph + ((1 + nc)*om*t0)/nc)))/(4.*(-1 + nc**2)*om**2*dsqrt(1 + xi**2))
    end if

    return;
end function pulse_alpha_z_sin2




double complex function pulse_alpha_x_sin2( t_ ) result(alpha_x)
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

    if( dabs( nc - 1d0 ) < 1d-6 ) then

        !        alpha_x = -(E0*xi*(cdcos(2*om*t) - cdcos(2*om*t0) + 4*cdsin(om*t) &
        !              - 4*cdsin(om*t0)))/(8.*om**2*dsqrt(1 + xi**2))
        
        alpha_x =   -(E0*xi*(cdcos(ph + 2*om*t) - cdcos(ph + 2*om*t0) + 2*om*t*dsin(ph) &
              - 2*om*t0*dsin(ph) + 4*cdsin(ph + om*t) - 4*cdsin(ph + om*t0)))/(8.*om**2*dsqrt(1 + xi**2))
        
    else

        !        alpha_x = -(E0*xi*(-(nc*(1 + nc)*cdcos(((-1 + nc)*(Pi - 2*om*t))/(2.*nc))) + &
        !              nc*(1 + nc)*cdcos(((-1 + nc)*(Pi - 2*om*t0))/(2.*nc)) + &
        !              (-1 + nc)*(2*(1 + nc)*cdsin(om*t) - 2*(1 + nc)*cdsin(om*t0) + &
        !              nc*(cdsin((Pi - 2*(1 + nc)*om*t)/(2.*nc)) - &
        !              cdsin((Pi - 2*(1 + nc)*om*t0)/(2.*nc))))))/(4.* &
        !              (-1 + nc**2)*om**2*dsqrt(1 + xi**2))
        alpha_x = -(E0*xi*(2*(-1 + nc**2)*cdsin(ph + om*t) - nc*(1 + nc)*cdsin(ph &
              + (Pi + 2*(-1 + nc)*om*t)/(2.*nc)) + (-1 + nc)*nc*cdsin((Pi - 2*(nc*ph &
              + (1 + nc)*om*t))/(2.*nc)) - 2*(-1 + nc**2)*cdsin(ph + om*t0) + nc*(1 + &
              nc)*cdsin(ph + (Pi + 2*(-1 + nc)*om*t0)/(2.*nc)) + (-1 + nc)*nc*cdsin(ph - &
              (Pi - 2*(1 + nc)*om*t0)/(2.*nc))))/(4.*(-1 + nc**2)*om**2*dsqrt(1 + xi**2));
    end if

    return;
end function pulse_alpha_x_sin2
