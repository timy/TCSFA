double complex function pulse_A_z_sin2( t ) result(f)
    use mod_pulse;
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = -((E0*cdsin(om*t)*cdsin((om*t)/(2.*nc))**2)/om)
    else if( t_ > Tp) then
        f = dcmplx(0d0)
    end if
    return;
end function pulse_A_z_sin2




double complex function pulse_A_x_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_A_x_sin2




double complex function pulse_E_z_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
!        f = (E0*(nc*cdcos(om*t)*(-1 + cdcos((om*t)/nc)) + cdsin(om*t)*cdsin((om*t)/nc)))/(2.*nc)
        f = (E0*cdsin((om*t)/(2.*nc))*(cdcos((om*t)/(2.*nc))*cdsin(om*t) + nc*cdcos(om*t)*cdsin((om*t)/(2.*nc))))/nc;
    else if( t_ > Tp) then
        f = dcmplx(0d0)
    end if


    return;
end function pulse_E_z_sin2




double complex function pulse_E_x_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_E_x_sin2




double complex function pulse_alpha_z_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = -(E0*(-1 + cdcos(om*t)*(1 - nc**2 + nc**2*cdcos((om*t)/nc)) + nc*cdsin(om*t)*cdsin((om*t)/nc)))/(2.*(-1 + nc**2)*om**2)
    else if( t_ > Tp) then
        f = dcmplx(-(E0*(-1 + dcos(om*Tp)*(1 - nc**2 + nc**2*dcos((om*Tp)/nc)) + nc*dsin(om*Tp)*dsin((om*Tp)/nc)))/(2.*(-1 + nc**2)*om**2))
    end if


    return;
end function pulse_alpha_z_sin2




double complex function pulse_alpha_x_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_alpha_x_sin2



double complex function pulse_alpha2_z_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = (E0**2*(12*om*t - 6*cdsin(2*om*t) + nc*((8*cdsin((2 + 1/nc)*om*t))/(1 + 2*nc) - 16*cdsin((om*t)/nc) + 2*cdsin((2*om*t)/nc) - cdsin((2*(-1 + nc)*om*t)/nc)/(-1 + nc) - cdsin((2*(1 + nc)*om*t)/nc)/(1 + nc) + (8*cdsin(((-1 + 2*nc)*om*t)/nc))/(-1 + 2*nc))))/(64.*om**3)
    else if( t_ > Tp) then
        f = dcmplx((E0**2*(12*om*Tp - 6*dsin(2*om*Tp) + nc*((8*dsin((2 + 1/nc)*om*Tp))/(1 + 2*nc) - 16*dsin((om*Tp)/nc) + 2*dsin((2*om*Tp)/nc) - dsin((2*(-1 + nc)*om*Tp)/nc)/(-1 + nc) - dsin((2*(1 + nc)*om*Tp)/nc)/(1 + nc) + (8*dsin(((-1 + 2*nc)*om*Tp)/nc))/(-1 + 2*nc))))/(64.*om**3))
    end if


    return;
end function pulse_alpha2_z_sin2




double complex function pulse_alpha2_x_sin2( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    f = dcmplx(0d0)
    return;
end function pulse_alpha2_x_sin2
