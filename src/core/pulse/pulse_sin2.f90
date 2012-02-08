#include '../../include/inc_field.h'
! ////////////////////////////////////////////////////////////////////////////////
! the vector potential in z-direction
double complex function pulse_A_z_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = -((E0*cdsin(OM*t)*cdsin((OM*t)/(2.*NC))**2)/OM)
    else if( t_ > Tp) then
        f = dcmplx(0d0)
    end if
    return;
end function pulse_A_z_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the vector potential in x-direction
double complex function pulse_A_x_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_A_x_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the electric field in z-direction
double complex function pulse_E_z_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
!        f = (E0*(nc*cdcos(om*t)*(-1 + cdcos((om*t)/nc)) + cdsin(om*t)*cdsin((om*t)/nc)))/(2.*nc)
        f = (E0*cdsin((OM*t)/(2.*NC))*(cdcos((OM*t)/(2.*NC))*cdsin(OM*t) + NC*cdcos(OM*t)*cdsin((OM*t)/(2.*NC))))/NC;
    else if( t_ > Tp) then
        f = dcmplx(0d0)
    end if
    return;
end function pulse_E_z_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the electric field in x-direction
double complex function pulse_E_x_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_E_x_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the alpha in z-direction
double complex function pulse_alpha_z_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = -(E0*(-1 + cdcos(OM*t)*(1 - NC**2 + NC**2*cdcos((OM*t)/NC)) + NC*cdsin(OM*t)*cdsin((OM*t)/NC)))/(2.*(-1 + NC**2)*OM**2)
    else if( t_ > Tp) then
        f = dcmplx(-(E0*(-1 + dcos(OM*Tp)*(1 - NC**2 + NC**2*dcos((OM*Tp)/NC)) + NC*dsin(OM*Tp)*dsin((OM*Tp)/NC)))/(2.*(-1 + NC**2)*OM**2))
    end if
    return;
end function pulse_alpha_z_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the alpha in x-direction
double complex function pulse_alpha_x_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    f = dcmplx(0d0)
    return;
end function pulse_alpha_x_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the alpha^2 in z-direction
double complex function pulse_alpha2_z_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    t_ = dreal(t)
    if( t_ <= 0d0 ) then
        f = dcmplx(0d0)
    else if( t_ > 0d0 .and. t_ < Tp ) then
        f = (E0**2*(12*OM*t - 6*cdsin(2*OM*t) + NC*((8*cdsin((2 + 1/NC)*OM*t))/(1 + 2*NC) - 16*cdsin((OM*t)/NC) + 2*cdsin((2*OM*t)/NC) - cdsin((2*(-1 + NC)*OM*t)/NC)/(-1 + NC) - cdsin((2*(1 + NC)*OM*t)/NC)/(1 + NC) + (8*cdsin(((-1 + 2*NC)*OM*t)/NC))/(-1 + 2*NC))))/(64.*OM**3)
    else if( t_ > Tp) then
        f = dcmplx((E0**2*(12*OM*Tp - 6*dsin(2*OM*Tp) + NC*((8*dsin((2 + 1/NC)*OM*Tp))/(1 + 2*NC) - 16*dsin((OM*Tp)/NC) + 2*dsin((2*OM*Tp)/NC) - dsin((2*(-1 + NC)*OM*Tp)/NC)/(-1 + NC) - dsin((2*(1 + NC)*OM*Tp)/NC)/(1 + NC) + (8*dsin(((-1 + 2*NC)*OM*Tp)/NC))/(-1 + 2*NC))))/(64.*OM**3))
    end if


    return;
end function pulse_alpha2_z_sin2

! ////////////////////////////////////////////////////////////////////////////////
! the alpha^2 in x-direction
double complex function pulse_alpha2_x_sin2( t ) result(f)
    implicit none;
    double complex, intent(in):: t;
    double precision:: t_
    f = dcmplx(0d0)
    return;
end function pulse_alpha2_x_sin2
