#include '../include/inc_atom.h'

double complex function action_W_im( ts )
    use mod_p0
    use mod_pulse, only: om, E0, nc, xi, Pi

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
    action_W_im = Ip*t0 - Ip*ts + (64*om**3*p0x**2*t0 - 320*nc**2*om**3*p0x**2*t0 &
         + 256*nc**4*om**3*p0x**2*t0 - 64*om**3*p0x**2*ts + 320*nc**2*om**3*p0x**2*ts &
         - 256*nc**4*om**3*p0x**2*ts + 12*E0**2*om*t0*xi**2 - 60*E0**2*nc**2*om*t0*xi**2 &
         + 48*E0**2*nc**4*om*t0*xi**2 + 64*om**3*p0x**2*t0*xi**2 - 320*nc**2*om**3*p0x**2*t0*xi**2 &
         + 256*nc**4*om**3*p0x**2*t0*xi**2 - 12*E0**2*om*ts*xi**2 + 60*E0**2*nc**2*om*ts*xi**2 &
         - 48*E0**2*nc**4*om*ts*xi**2 - 64*om**3*p0x**2*ts*xi**2 + 320*nc**2*om**3*p0x**2*ts*xi**2 &
         - 256*nc**4*om**3*p0x**2*ts*xi**2 + 32*E0*nc*(-1 - nc + 4*nc**2 + 4*nc**3)*om*p0x*xi*dsqrt(1 &
         + xi**2)*cdcos(((-1 + nc)*(Pi - 2*om*t0))/(2.*nc)) - 32*E0*nc*(-1 - nc + 4*nc**2 &
         + 4*nc**3)*om*p0x*xi*dsqrt(1 + xi**2)*cdcos(((-1 + nc)*(Pi - 2*om*ts))/(2.*nc)) - &
         64*E0*(1 - 5*nc**2 + 4*nc**4)*om*p0x*xi*dsqrt(1 + xi**2)*cdsin(om*t0) + &
         6*E0**2*xi**2*cdsin(2*om*t0) - 30*E0**2*nc**2*xi**2*cdsin(2*om*t0) + &
         24*E0**2*nc**4*xi**2*cdsin(2*om*t0) + 16*E0**2*nc*xi**2*cdsin((Pi - 2*om*t0)/(2.*nc)) &
         - 80*E0**2*nc**3*xi**2*cdsin((Pi - 2*om*t0)/(2.*nc)) + 64*E0**2*nc**5*xi**2*cdsin((Pi &
         - 2*om*t0)/(2.*nc)) - 2*E0**2*nc*xi**2*cdsin((Pi - 2*om*t0)/nc) + &
         10*E0**2*nc**3*xi**2*cdsin((Pi - 2*om*t0)/nc) - 8*E0**2*nc**5*xi**2*cdsin((Pi - 2*om*t0)/nc) &
         - E0**2*nc*xi**2*cdsin((Pi + 2*(-1 + nc)*om*t0)/nc) - E0**2*nc**2*xi**2*cdsin((Pi + 2*(-1 &
         + nc)*om*t0)/nc) +  4*E0**2*nc**3*xi**2*cdsin((Pi + 2*(-1 + nc)*om*t0)/nc) + &
         4*E0**2*nc**4*xi**2*cdsin((Pi + 2*(-1 + nc)*om*t0)/nc) + 8*E0**2*nc*xi**2*cdsin((Pi - 2*om*t0 &
         + 4*nc*om*t0)/(2.*nc)) + 16*E0**2*nc**2*xi**2*cdsin((Pi - 2*om*t0 + 4*nc*om*t0)/(2.*nc)) &
         - 8*E0**2*nc**3*xi**2*cdsin((Pi - 2*om*t0 + 4*nc*om*t0)/(2.*nc)) - &
         16*E0**2*nc**4*xi**2*cdsin((Pi - 2*om*t0 + 4*nc*om*t0)/(2.*nc)) - 32*E0*nc*om*p0x*xi*dsqrt(1 &
         + xi**2)*cdsin((Pi - 2*(1 + nc)*om*t0)/(2.*nc)) + 32*E0*nc**2*om*p0x*xi*dsqrt(1 + &
         xi**2)*cdsin((Pi - 2*(1 + nc)*om*t0)/(2.*nc)) + 128*E0*nc**3*om*p0x*xi*dsqrt(1 + &
         xi**2)*cdsin((Pi - 2*(1 + nc)*om*t0)/(2.*nc)) - 128*E0*nc**4*om*p0x*xi*dsqrt(1 + &
         xi**2)*cdsin((Pi - 2*(1 + nc)*om*t0)/(2.*nc)) - E0**2*nc*xi**2*cdsin((Pi - 2*(1 + &
         nc)*om*t0)/nc) + E0**2*nc**2*xi**2*cdsin((Pi - 2*(1 + nc)*om*t0)/nc) + &
         4*E0**2*nc**3*xi**2*cdsin((Pi - 2*(1 + nc)*om*t0)/nc) - 4*E0**2*nc**4*xi**2*cdsin((Pi - &
         2*(1 + nc)*om*t0)/nc) + 8*E0**2*nc*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*t0)/(2.*nc)) - &
         16*E0**2*nc**2*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*t0)/(2.*nc)) - &
         8*E0**2*nc**3*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*t0)/(2.*nc)) + &
         16*E0**2*nc**4*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*t0)/(2.*nc)) + 64*E0*(1 - 5*nc**2 + &
         4*nc**4)*om*p0x*xi*dsqrt(1 + xi**2)*cdsin(om*ts) - 6*E0**2*xi**2*cdsin(2*om*ts) + &
         30*E0**2*nc**2*xi**2*cdsin(2*om*ts) - 24*E0**2*nc**4*xi**2*cdsin(2*om*ts) + &
         2*(6*E0**2*om*t0 - 30*E0**2*nc**2*om*t0 + 24*E0**2*nc**4*om*t0 + 32*om**3*p0z**2*t0 - &
         160*nc**2*om**3*p0z**2*t0 + 128*nc**4*om**3*p0z**2*t0 - 6*E0**2*om*ts + 30*E0**2*nc**2*om*ts &
         - 24*E0**2*nc**4*om*ts - 32*om**3*p0z**2*ts + 160*nc**2*om**3*p0z**2*ts - &
         128*nc**4*om**3*p0z**2*ts + 32*om**3*p0z**2*t0*xi**2 - 160*nc**2*om**3*p0z**2*t0*xi**2 + &
         128*nc**4*om**3*p0z**2*t0*xi**2 - 32*om**3*p0z**2*ts*xi**2 + 160*nc**2*om**3*p0z**2*ts*xi**2 - &
         128*nc**4*om**3*p0z**2*ts*xi**2 - 32*E0*(-1 + 4*nc**2)*om*p0z*dsqrt(1 + xi**2)*cdcos(om*t0)*(1 - &
         nc**2 + nc**2*cdcos((om*t0)/nc)) + 32*E0*(-1 + 4*nc**2)*om*p0z*dsqrt(1 + xi**2)*cdcos(om*ts)*(1 - &
         nc**2 + nc**2*cdcos((om*ts)/nc)) + E0**2*(16*nc**2*(-1 + nc**2)*cdcos((om*t0)/nc) - &
         (-1 + 4*nc**2)*(-3 + 3*nc**2 + nc**2*cdcos((2*om*t0)/nc)))*cdsin(2*om*t0) - &
         8*E0**2*nc*cdsin((om*t0)/nc) + 40*E0**2*nc**3*cdsin((om*t0)/nc) - &
         32*E0**2*nc**5*cdsin((om*t0)/nc) + 2*E0**2*nc*cdcos(om*t0)**2*(4 - 4*nc**2 + (-1 + &
         4*nc**2)*cdcos((om*t0)/nc))*cdsin((om*t0)/nc) + 32*E0*nc*om*p0z*dsqrt(1 + &
         xi**2)*cdsin(om*t0)*cdsin((om*t0)/nc) - 128*E0*nc**3*om*p0z*dsqrt(1 + &
         xi**2)*cdsin(om*t0)*cdsin((om*t0)/nc) - 8*E0**2*nc*cdsin(om*t0)**2*cdsin((om*t0)/nc) + &
         8*E0**2*nc**3*cdsin(om*t0)**2*cdsin((om*t0)/nc) + E0**2*nc*cdsin((2*om*t0)/nc) - &
         5*E0**2*nc**3*cdsin((2*om*t0)/nc) + 4*E0**2*nc**5*cdsin((2*om*t0)/nc) + &
         E0**2*nc*cdsin(om*t0)**2*cdsin((2*om*t0)/nc) - 4*E0**2*nc**3*cdsin(om*t0)**2*cdsin((2*om*t0)/nc) &
         - E0**2*(16*nc**2*(-1 + nc**2)*cdcos((om*ts)/nc) - (-1 + 4*nc**2)*(-3 + 3*nc**2 + &
         nc**2*cdcos((2*om*ts)/nc)))*cdsin(2*om*ts) + 8*E0**2*nc*cdsin((om*ts)/nc) - &
         40*E0**2*nc**3*cdsin((om*ts)/nc) + 32*E0**2*nc**5*cdsin((om*ts)/nc) - &
         2*E0**2*nc*cdcos(om*ts)**2*(4 - 4*nc**2 + (-1 + 4*nc**2)*cdcos((om*ts)/nc))*cdsin((om*ts)/nc) - &
         32*E0*nc*om*p0z*dsqrt(1 + xi**2)*cdsin(om*ts)*cdsin((om*ts)/nc) + 128*E0*nc**3*om*p0z*dsqrt(1 &
         + xi**2)*cdsin(om*ts)*cdsin((om*ts)/nc) + 8*E0**2*nc*cdsin(om*ts)**2*cdsin((om*ts)/nc) &
         - 8*E0**2*nc**3*cdsin(om*ts)**2*cdsin((om*ts)/nc) - E0**2*nc*cdsin((2*om*ts)/nc) + &
         5*E0**2*nc**3*cdsin((2*om*ts)/nc) - 4*E0**2*nc**5*cdsin((2*om*ts)/nc) - &
         E0**2*nc*cdsin(om*ts)**2*cdsin((2*om*ts)/nc) + 4*E0**2*nc**3*cdsin(om*ts)**2*cdsin((2*om*ts)/nc)) &
         - 16*E0**2*nc*xi**2*cdsin((Pi - 2*om*ts)/(2.*nc)) + 80*E0**2*nc**3*xi**2*cdsin((Pi - &
         2*om*ts)/(2.*nc)) - 64*E0**2*nc**5*xi**2*cdsin((Pi - 2*om*ts)/(2.*nc)) + &
         2*E0**2*nc*xi**2*cdsin((Pi - 2*om*ts)/nc) - 10*E0**2*nc**3*xi**2*cdsin((Pi - 2*om*ts)/nc) + &
         8*E0**2*nc**5*xi**2*cdsin((Pi - 2*om*ts)/nc) + E0**2*nc*xi**2*cdsin((Pi + 2*(-1 + nc)*om*ts)/nc) &
         + E0**2*nc**2*xi**2*cdsin((Pi + 2*(-1 + nc)*om*ts)/nc) - 4*E0**2*nc**3*xi**2*cdsin((Pi + 2*(-1 &
         + nc)*om*ts)/nc) - 4*E0**2*nc**4*xi**2*cdsin((Pi + 2*(-1 + nc)*om*ts)/nc) - &
         8*E0**2*nc*xi**2*cdsin((Pi - 2*om*ts + 4*nc*om*ts)/(2.*nc)) - 16*E0**2*nc**2*xi**2*cdsin((Pi - &
         2*om*ts + 4*nc*om*ts)/(2.*nc)) + 8*E0**2*nc**3*xi**2*cdsin((Pi - 2*om*ts + 4*nc*om*ts)/(2.*nc)) &
         + 16*E0**2*nc**4*xi**2*cdsin((Pi - 2*om*ts + 4*nc*om*ts)/(2.*nc)) + 32*E0*nc*om*p0x*xi*dsqrt(1 + &
         xi**2)*cdsin((Pi - 2*(1 + nc)*om*ts)/(2.*nc)) - 32*E0*nc**2*om*p0x*xi*dsqrt(1 + xi**2)*cdsin((Pi &
         - 2*(1 + nc)*om*ts)/(2.*nc)) - 128*E0*nc**3*om*p0x*xi*dsqrt(1 + xi**2)*cdsin((Pi - 2*(1 + &
         nc)*om*ts)/(2.*nc)) + 128*E0*nc**4*om*p0x*xi*dsqrt(1 + xi**2)*cdsin((Pi - 2*(1 + &
         nc)*om*ts)/(2.*nc)) + E0**2*nc*xi**2*cdsin((Pi - 2*(1 + nc)*om*ts)/nc) - &
         E0**2*nc**2*xi**2*cdsin((Pi - 2*(1 + nc)*om*ts)/nc) - 4*E0**2*nc**3*xi**2*cdsin((Pi - &
         2*(1 + nc)*om*ts)/nc) + 4*E0**2*nc**4*xi**2*cdsin((Pi - 2*(1 + nc)*om*ts)/nc) - &
         8*E0**2*nc*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*ts)/(2.*nc)) + 16*E0**2*nc**2*xi**2*cdsin((Pi &
         - 2*(1 + 2*nc)*om*ts)/(2.*nc)) + 8*E0**2*nc**3*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*ts)/(2.*nc)) &
         - 16*E0**2*nc**4*xi**2*cdsin((Pi - 2*(1 + 2*nc)*om*ts)/(2.*nc)))/(128.*(1 - 5*nc**2 + &
         4*nc**4)*om**3*(1 + xi**2))

!!$     ! this one is for constant envelope
!!$     action_W_im = (-8*Ip*t0 - 4*p0x**2*t0 - 4*p0z**2*t0 + 8*Ip*ts + 4*p0x**2*ts + 4*p0z**2*ts + &
!!$           (8*E0*p0z*(cdsin(om*t0) - cdsin(om*ts)))/(om**2*dsqrt(1 + xi**2)) + (8*E0*p0x*xi*(cdsin(om*t0) &
!!$           - cdsin(om*ts)))/(om**2*dsqrt(1 + xi**2)) - (E0**2*(2*om*(t0 - ts) + cdsin(2*om*t0) - cdsin(2*om*ts)))/om**3)/8.
!!$     action_W_im = - action_W_im

!!$    action_w_im = ( 2d0 * ( t0 - ts ) * om * ( E0*E0 + 2d0 * ( 2d0 * Ip &
!!$          + p0x*p0x + p0z*p0z) * om * om ) + E0 * ( 8d0 * p0z * om &
!!$          *( cdcos( om * t0 ) - cdcos( om * ts ) ) + E0 &
!!$          * ( -cdsin( 2d0*om*t0 ) + cdsin( 2d0*om*ts ) ) ) ) &
!!$          / ( 8d0 * om * om * om )

    return;
end function action_W_im


