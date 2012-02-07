double complex function envelope_f( t ) result(f)
    use mod_pulse
    implicit none
    double complex, intent(in):: t
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ <= head ) then
        f = dcmplx(t/head)
    else if( t_ > head .and. t_ <= tail ) then
        f = dcmplx(1d0)
    else if( t_ > tail .and. t_ <= Tp ) then
        f = (Tp - t) / (Tp - tail)
    else ! t > Tp
        f = dcmplx(0d0)
    end if

    return
end function envelope_f


double complex function envelope_df( t ) result(f)
    use mod_pulse
    implicit none
    double complex, intent(in):: t
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ <= head ) then
        f = dcmplx(1d0/head)
    else if( t_ > head .and. t_ <= tail ) then
        f = dcmplx(0d0)
    else if( t_ > tail .and. t_ <= Tp ) then
        f = dcmplx(1d0/(tail-Tp))
    else ! t > Tp
        f = dcmplx(0d0)
    end if

    return
end function envelope_df

! vector potential A(t)
double complex function pulse_A_z_ramp( t ) result(A_z)
    use mod_pulse;
    implicit none;
    double complex, external:: envelope_f
    double complex, intent(in):: t;
    A_z = dcmplx( - E0 * envelope_f(t) * cdsin(om*t) / om )
    return;
end function pulse_A_z_ramp

double complex function pulse_A_x_ramp( t ) result(A_x)
    use mod_pulse;
    implicit none;
    double complex, intent(in):: t;
    A_x = dcmplx(0d0)
    return;    
end function pulse_A_x_ramp


! electric field E(t)
double complex function pulse_E_z_ramp( t ) result(E_z)
    use mod_pulse;
    implicit none;
    double complex, intent(in):: t;
    double complex, external:: envelope_f, envelope_df
    E_z = E0 * cdcos(om*t) * envelope_f(t) + E0 * cdsin(om*t) * envelope_df(t) / om
    return;
end function pulse_E_z_ramp

double complex function pulse_E_x_ramp( t ) result(E_x)
    use mod_pulse;
    implicit none;
    double complex, intent(in):: t;
    E_x = dcmplx(0d0)
    return;
end function pulse_E_x_ramp


! alpha
double complex function pulse_alpha_z_ramp( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex:: t;
    double precision:: t_

    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ <= head ) then
        f = (E0*(om*t*cdcos(om*t) - cdsin(om*t)))/(head*om**3)
    else if( t_ > head .and. t_ <= tail ) then
        f = (E0*(head*om*cdcos(om*t) - dsin(head*om)))/(head*om**3)
    else if( t_ > tail .and. t_ <= Tp ) then
        f = (E0*(head*om*(t - Tp)*cdcos(om*t) + (-tail + Tp)*dsin(head*om) + head*(-cdsin(om*t) + dsin(om*tail))))/(head*om**3*(tail - Tp))
    else ! t > Tp
        f = (E0*((Tp-tail)*dsin(head*om) + head*(dsin(om*tail) - dsin(om*Tp))))/(head*om**3*(tail - Tp))
    end if

    return;
end function pulse_alpha_z_ramp

double complex function pulse_alpha_x_ramp( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_alpha_x_ramp

! alpha2
double complex function pulse_alpha2_z_ramp( t ) result(f)
    use mod_pulse;

    implicit none;
    double complex:: t;
    double precision:: t_

    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ <= head ) then
        f = (E0**2*(4*(om*t)**3 + 3*cdsin(2*om*t) - 6*om*t*(cdcos(2*om*t) + om*t*cdsin(2*om*t))))/(24.*head**2*om**5)
    else if( t_ > head .and. t_ <= tail ) then
        f = (E0**2*(4*head**3*om**3 + 3*dsin(2*head*om) - 6*head*om*(dcos(2*head*om) + head*om*dsin(2*head*om))))/(24.*head**2*om**5) + (E0**2*(2*om*(-head + t) + dsin(2*head*om) - cdsin(2*om*t)))/(4.*om**3)
    else if( t_ > tail .and. t_ <= Tp ) then
        f = (E0**2*((4*head**3*om**3 + 3*dsin(2*head*om) - 6*head*om*(dcos(2*head*om) + head*om*dsin(2*head*om)))/head**2 + 6*om**2*(2*om*(-head + tail) + dsin(2*head*om) - dsin(2*om*tail)) + (4*om**3*t*(t**2 - 3*t*Tp + 3*Tp**2) - 4*om**3*tail*(tail**2 - 3*tail*Tp + 3*Tp**2) +  6*om*(-t + Tp)*cdcos(2*om*t) + 6*om*(tail - Tp)*dcos(2*om*tail) + (3 - 6*om**2*(t - Tp)**2)*cdsin(2*om*t) - (3 - 6*om**2*(tail - Tp)**2)*dsin(2*om*tail))/(tail - Tp)**2))/(24.*om**5)
    else ! t > Tp
        f = (E0**2*((4*head**3*om**3 + 3*dsin(2*head*om) - 6*head*om*(dcos(2*head*om) + head*om*dsin(2*head*om)))/head**2 + 6*om**2*(2*om*(-head + tail) + dsin(2*head*om) - dsin(2*om*tail)) + (-4*om**3*(tail - Tp)**3 - 3*dsin(2*om*tail) + 6*om*(tail - Tp)*(dcos(2*om*tail) + om*(tail - Tp)*dsin(2*om*tail)) + 3*dsin(2*om*Tp))/(tail - Tp)**2))/(24.*om**5)
    end if

    return;
end function pulse_alpha2_z_ramp

double complex function pulse_alpha2_x_ramp( t ) result(f)
    use mod_pulse;
    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_alpha2_x_ramp
