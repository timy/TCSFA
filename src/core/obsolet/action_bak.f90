#include '../include/inc_atom.h'
#include '../include/inc_field.h'
#include '../include/inc_const.h'

double complex function action_W_im( ts )
    use mod_p0
    implicit none;
    double complex, intent(in):: ts;
    double complex, external:: v2_integrand, cedint;
    double complex:: t0;
    double precision, parameter:: Ip = IONIZATION_IP
    double precision:: p0x, p0z
#ifdef IM_PLOT_INTEGRAND
    integer:: i
    double complex:: t
#endif

    t0 = dcmplx( dreal(ts), 0d0 );
    p0x = p0_x;
    p0z = p0_z;
    if ( dimag(ts)+1d-6 <= 0d0 ) stop 'error: Im[t] < 0';


#ifdef IM_PLOT_INTEGRAND
    open(IM_PLOT_FILE_ID, file=IM_PLOT_FILE_NAME)
    do i = 1, IM_PLOT_N_PTS
        t = dcmplx( dreal(ts), ( dimag(ts) + IM_PLOT_OFFSET ) / (IM_PLOT_N_PTS-1) * (i - 1) );
        write(IM_PLOT_FILE_ID, '(3(e15.8,2x))'), dimag(t), v2_integrand( t ) + Ip;
    end do
    close(IM_PLOT_FILE_ID)
#endif ! IM_PLOT_INTEGRAND

    ! pay attention here that the px has been set as 0 for dieter's test (sub-barrier correction)
    ! p0x = 0d0
    ! end attention

    ! this one is for sin^2
    action_W_im = Ip*t0 - Ip*ts + (64*OM**3*p0x**2*t0 - 320*NC**2*OM**3*p0x**2*t0 &
         + 256*NC**4*OM**3*p0x**2*t0 - 64*OM**3*p0x**2*ts + 320*NC**2*OM**3*p0x**2*ts &
         - 256*NC**4*OM**3*p0x**2*ts + 12*E0**2*OM*t0*XI**2 - 60*E0**2*NC**2*OM*t0*XI**2 &
         + 48*E0**2*NC**4*OM*t0*XI**2 + 64*OM**3*p0x**2*t0*XI**2 - 320*NC**2*OM**3*p0x**2*t0*XI**2 &
         + 256*NC**4*OM**3*p0x**2*t0*XI**2 - 12*E0**2*OM*ts*XI**2 + 60*E0**2*NC**2*OM*ts*XI**2 &
         - 48*E0**2*NC**4*OM*ts*XI**2 - 64*OM**3*p0x**2*ts*XI**2 + 320*NC**2*OM**3*p0x**2*ts*XI**2 &
         - 256*NC**4*OM**3*p0x**2*ts*XI**2 + 32*E0*NC*(-1 - NC + 4*NC**2 + 4*NC**3)*OM*p0x*XI*dsqrt(1 &
         + XI**2)*cdcos(((-1 + NC)*(PI - 2*OM*t0))/(2.*NC)) - 32*E0*NC*(-1 - NC + 4*NC**2 &
         + 4*NC**3)*OM*p0x*XI*dsqrt(1 + XI**2)*cdcos(((-1 + NC)*(PI - 2*OM*ts))/(2.*NC)) - &
         64*E0*(1 - 5*NC**2 + 4*NC**4)*OM*p0x*XI*dsqrt(1 + XI**2)*cdsin(OM*t0) + &
         6*E0**2*XI**2*cdsin(2*OM*t0) - 30*E0**2*NC**2*XI**2*cdsin(2*OM*t0) + &
         24*E0**2*NC**4*XI**2*cdsin(2*OM*t0) + 16*E0**2*NC*XI**2*cdsin((PI - 2*OM*t0)/(2.*NC)) &
         - 80*E0**2*NC**3*XI**2*cdsin((PI - 2*OM*t0)/(2.*NC)) + 64*E0**2*NC**5*XI**2*cdsin((PI &
         - 2*OM*t0)/(2.*NC)) - 2*E0**2*NC*XI**2*cdsin((PI - 2*OM*t0)/NC) + &
         10*E0**2*NC**3*XI**2*cdsin((PI - 2*OM*t0)/NC) - 8*E0**2*NC**5*XI**2*cdsin((PI - 2*OM*t0)/NC) &
         - E0**2*NC*XI**2*cdsin((PI + 2*(-1 + NC)*OM*t0)/NC) - E0**2*NC**2*XI**2*cdsin((PI + 2*(-1 &
         + NC)*OM*t0)/NC) +  4*E0**2*NC**3*XI**2*cdsin((PI + 2*(-1 + NC)*OM*t0)/NC) + &
         4*E0**2*NC**4*XI**2*cdsin((PI + 2*(-1 + NC)*OM*t0)/NC) + 8*E0**2*NC*XI**2*cdsin((PI - 2*OM*t0 &
         + 4*NC*OM*t0)/(2.*NC)) + 16*E0**2*NC**2*XI**2*cdsin((PI - 2*OM*t0 + 4*NC*OM*t0)/(2.*NC)) &
         - 8*E0**2*NC**3*XI**2*cdsin((PI - 2*OM*t0 + 4*NC*OM*t0)/(2.*NC)) - &
         16*E0**2*NC**4*XI**2*cdsin((PI - 2*OM*t0 + 4*NC*OM*t0)/(2.*NC)) - 32*E0*NC*OM*p0x*XI*dsqrt(1 &
         + XI**2)*cdsin((PI - 2*(1 + NC)*OM*t0)/(2.*NC)) + 32*E0*NC**2*OM*p0x*XI*dsqrt(1 + &
         XI**2)*cdsin((PI - 2*(1 + NC)*OM*t0)/(2.*NC)) + 128*E0*NC**3*OM*p0x*XI*dsqrt(1 + &
         XI**2)*cdsin((PI - 2*(1 + NC)*OM*t0)/(2.*NC)) - 128*E0*NC**4*OM*p0x*XI*dsqrt(1 + &
         XI**2)*cdsin((PI - 2*(1 + NC)*OM*t0)/(2.*NC)) - E0**2*NC*XI**2*cdsin((PI - 2*(1 + &
         NC)*OM*t0)/NC) + E0**2*NC**2*XI**2*cdsin((PI - 2*(1 + NC)*OM*t0)/NC) + &
         4*E0**2*NC**3*XI**2*cdsin((PI - 2*(1 + NC)*OM*t0)/NC) - 4*E0**2*NC**4*XI**2*cdsin((PI - &
         2*(1 + NC)*OM*t0)/NC) + 8*E0**2*NC*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*t0)/(2.*NC)) - &
         16*E0**2*NC**2*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*t0)/(2.*NC)) - &
         8*E0**2*NC**3*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*t0)/(2.*NC)) + &
         16*E0**2*NC**4*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*t0)/(2.*NC)) + 64*E0*(1 - 5*NC**2 + &
         4*NC**4)*OM*p0x*XI*dsqrt(1 + XI**2)*cdsin(OM*ts) - 6*E0**2*XI**2*cdsin(2*OM*ts) + &
         30*E0**2*NC**2*XI**2*cdsin(2*OM*ts) - 24*E0**2*NC**4*XI**2*cdsin(2*OM*ts) + &
         2*(6*E0**2*OM*t0 - 30*E0**2*NC**2*OM*t0 + 24*E0**2*NC**4*OM*t0 + 32*OM**3*p0z**2*t0 - &
         160*NC**2*OM**3*p0z**2*t0 + 128*NC**4*OM**3*p0z**2*t0 - 6*E0**2*OM*ts + 30*E0**2*NC**2*OM*ts &
         - 24*E0**2*NC**4*OM*ts - 32*OM**3*p0z**2*ts + 160*NC**2*OM**3*p0z**2*ts - &
         128*NC**4*OM**3*p0z**2*ts + 32*OM**3*p0z**2*t0*XI**2 - 160*NC**2*OM**3*p0z**2*t0*XI**2 + &
         128*NC**4*OM**3*p0z**2*t0*XI**2 - 32*OM**3*p0z**2*ts*XI**2 + 160*NC**2*OM**3*p0z**2*ts*XI**2 - &
         128*NC**4*OM**3*p0z**2*ts*XI**2 - 32*E0*(-1 + 4*NC**2)*OM*p0z*dsqrt(1 + XI**2)*cdcos(OM*t0)*(1 - &
         NC**2 + NC**2*cdcos((OM*t0)/NC)) + 32*E0*(-1 + 4*NC**2)*OM*p0z*dsqrt(1 + XI**2)*cdcos(OM*ts)*(1 - &
         NC**2 + NC**2*cdcos((OM*ts)/NC)) + E0**2*(16*NC**2*(-1 + NC**2)*cdcos((OM*t0)/NC) - &
         (-1 + 4*NC**2)*(-3 + 3*NC**2 + NC**2*cdcos((2*OM*t0)/NC)))*cdsin(2*OM*t0) - &
         8*E0**2*NC*cdsin((OM*t0)/NC) + 40*E0**2*NC**3*cdsin((OM*t0)/NC) - &
         32*E0**2*NC**5*cdsin((OM*t0)/NC) + 2*E0**2*NC*cdcos(OM*t0)**2*(4 - 4*NC**2 + (-1 + &
         4*NC**2)*cdcos((OM*t0)/NC))*cdsin((OM*t0)/NC) + 32*E0*NC*OM*p0z*dsqrt(1 + &
         XI**2)*cdsin(OM*t0)*cdsin((OM*t0)/NC) - 128*E0*NC**3*OM*p0z*dsqrt(1 + &
         XI**2)*cdsin(OM*t0)*cdsin((OM*t0)/NC) - 8*E0**2*NC*cdsin(OM*t0)**2*cdsin((OM*t0)/NC) + &
         8*E0**2*NC**3*cdsin(OM*t0)**2*cdsin((OM*t0)/NC) + E0**2*NC*cdsin((2*OM*t0)/NC) - &
         5*E0**2*NC**3*cdsin((2*OM*t0)/NC) + 4*E0**2*NC**5*cdsin((2*OM*t0)/NC) + &
         E0**2*NC*cdsin(OM*t0)**2*cdsin((2*OM*t0)/NC) - 4*E0**2*NC**3*cdsin(OM*t0)**2*cdsin((2*OM*t0)/NC) &
         - E0**2*(16*NC**2*(-1 + NC**2)*cdcos((OM*ts)/NC) - (-1 + 4*NC**2)*(-3 + 3*NC**2 + &
         NC**2*cdcos((2*OM*ts)/NC)))*cdsin(2*OM*ts) + 8*E0**2*NC*cdsin((OM*ts)/NC) - &
         40*E0**2*NC**3*cdsin((OM*ts)/NC) + 32*E0**2*NC**5*cdsin((OM*ts)/NC) - &
         2*E0**2*NC*cdcos(OM*ts)**2*(4 - 4*NC**2 + (-1 + 4*NC**2)*cdcos((OM*ts)/NC))*cdsin((OM*ts)/NC) - &
         32*E0*NC*OM*p0z*dsqrt(1 + XI**2)*cdsin(OM*ts)*cdsin((OM*ts)/NC) + 128*E0*NC**3*OM*p0z*dsqrt(1 &
         + XI**2)*cdsin(OM*ts)*cdsin((OM*ts)/NC) + 8*E0**2*NC*cdsin(OM*ts)**2*cdsin((OM*ts)/NC) &
         - 8*E0**2*NC**3*cdsin(OM*ts)**2*cdsin((OM*ts)/NC) - E0**2*NC*cdsin((2*OM*ts)/NC) + &
         5*E0**2*NC**3*cdsin((2*OM*ts)/NC) - 4*E0**2*NC**5*cdsin((2*OM*ts)/NC) - &
         E0**2*NC*cdsin(OM*ts)**2*cdsin((2*OM*ts)/NC) + 4*E0**2*NC**3*cdsin(OM*ts)**2*cdsin((2*OM*ts)/NC)) &
         - 16*E0**2*NC*XI**2*cdsin((PI - 2*OM*ts)/(2.*NC)) + 80*E0**2*NC**3*XI**2*cdsin((PI - &
         2*OM*ts)/(2.*NC)) - 64*E0**2*NC**5*XI**2*cdsin((PI - 2*OM*ts)/(2.*NC)) + &
         2*E0**2*NC*XI**2*cdsin((PI - 2*OM*ts)/NC) - 10*E0**2*NC**3*XI**2*cdsin((PI - 2*OM*ts)/NC) + &
         8*E0**2*NC**5*XI**2*cdsin((PI - 2*OM*ts)/NC) + E0**2*NC*XI**2*cdsin((PI + 2*(-1 + NC)*OM*ts)/NC) &
         + E0**2*NC**2*XI**2*cdsin((PI + 2*(-1 + NC)*OM*ts)/NC) - 4*E0**2*NC**3*XI**2*cdsin((PI + 2*(-1 &
         + NC)*OM*ts)/NC) - 4*E0**2*NC**4*XI**2*cdsin((PI + 2*(-1 + NC)*OM*ts)/NC) - &
         8*E0**2*NC*XI**2*cdsin((PI - 2*OM*ts + 4*NC*OM*ts)/(2.*NC)) - 16*E0**2*NC**2*XI**2*cdsin((PI - &
         2*OM*ts + 4*NC*OM*ts)/(2.*NC)) + 8*E0**2*NC**3*XI**2*cdsin((PI - 2*OM*ts + 4*NC*OM*ts)/(2.*NC)) &
         + 16*E0**2*NC**4*XI**2*cdsin((PI - 2*OM*ts + 4*NC*OM*ts)/(2.*NC)) + 32*E0*NC*OM*p0x*XI*dsqrt(1 + &
         XI**2)*cdsin((PI - 2*(1 + NC)*OM*ts)/(2.*NC)) - 32*E0*NC**2*OM*p0x*XI*dsqrt(1 + XI**2)*cdsin((PI &
         - 2*(1 + NC)*OM*ts)/(2.*NC)) - 128*E0*NC**3*OM*p0x*XI*dsqrt(1 + XI**2)*cdsin((PI - 2*(1 + &
         NC)*OM*ts)/(2.*NC)) + 128*E0*NC**4*OM*p0x*XI*dsqrt(1 + XI**2)*cdsin((PI - 2*(1 + &
         NC)*OM*ts)/(2.*NC)) + E0**2*NC*XI**2*cdsin((PI - 2*(1 + NC)*OM*ts)/NC) - &
         E0**2*NC**2*XI**2*cdsin((PI - 2*(1 + NC)*OM*ts)/NC) - 4*E0**2*NC**3*XI**2*cdsin((PI - &
         2*(1 + NC)*OM*ts)/NC) + 4*E0**2*NC**4*XI**2*cdsin((PI - 2*(1 + NC)*OM*ts)/NC) - &
         8*E0**2*NC*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*ts)/(2.*NC)) + 16*E0**2*NC**2*XI**2*cdsin((PI &
         - 2*(1 + 2*NC)*OM*ts)/(2.*NC)) + 8*E0**2*NC**3*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*ts)/(2.*NC)) &
         - 16*E0**2*NC**4*XI**2*cdsin((PI - 2*(1 + 2*NC)*OM*ts)/(2.*NC)))/(128.*(1 - 5*NC**2 + &
         4*NC**4)*OM**3*(1 + XI**2))

!!$     ! this one is for constant envelope
!!$     action_W_im = (-8*Ip*t0 - 4*p0x**2*t0 - 4*p0z**2*t0 + 8*Ip*ts + 4*p0x**2*ts + 4*p0z**2*ts + &
!!$           (8*E0*p0z*(cdsin(OM*t0) - cdsin(OM*ts)))/(OM**2*dsqrt(1 + XI**2)) + (8*E0*p0x*XI*(cdsin(OM*t0) &
!!$           - cdsin(OM*ts)))/(OM**2*dsqrt(1 + XI**2)) - (E0**2*(2*OM*(t0 - ts) + cdsin(2*OM*t0) - cdsin(2*OM*ts)))/OM**3)/8.
!!$     action_W_im = - action_W_im

!!$    action_w_im = ( 2d0 * ( t0 - ts ) * OM * ( E0*E0 + 2d0 * ( 2d0 * Ip &
!!$          + p0x*p0x + p0z*p0z) * OM * OM ) + E0 * ( 8d0 * p0z * OM &
!!$          *( cdcos( OM * t0 ) - cdcos( OM * ts ) ) + E0 &
!!$          * ( -cdsin( 2d0*OM*t0 ) + cdsin( 2d0*OM*ts ) ) ) ) &
!!$          / ( 8d0 * OM * OM * OM )

    return;
end function action_W_im


