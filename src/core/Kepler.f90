double precision function asymptotic_angle( x0, y0, vx0, vy0, energy, charge ) result(angle)
    implicit none;
    double precision, intent(in):: x0, y0, vx0, vy0, energy, charge;
    double precision:: L, B, C, r0, theta_v0, theta_r0, theta_diff, T1, T2;
    double precision, parameter:: PI = 2d0 * asin(1d0);
    
    ! direction of the initial velocity
    if( vx0 < 0 ) then
        theta_v0 = datan( vy0 / vx0 ) + PI;
    else
        theta_v0 = datan( vy0 / vx0 );
    end if

    ! direction of the initial position
    if( x0 < 0 ) then
        theta_r0 = datan( y0 / x0 ) + PI;
    else
        theta_r0 = datan( y0 / x0 );
    end if

 

    ! if there is no long range central force
    if( dabs(charge) < 1d-8 ) then
        angle = theta_v0;
        return;
    end if

    L = x0 * vy0 - vx0 * y0;

    B = charge / ( L * L );
    C = dsqrt( B * B + 2d0 * energy / ( L * L ) );
    r0 = dsqrt( x0 * x0 + y0 * y0 );

    theta_diff = theta_v0 - theta_r0;
    if( theta_diff < 0 ) theta_diff = theta_diff + 2d0 * PI;


    T1 = dacos( -B / C );
    T2 = dacos( ( 1d0 - B * r0 ) / ( C * r0 ) );
    if( theta_diff >= 0d0 .and. theta_diff < 0.5d0*PI ) then
        angle = theta_r0 + T1 - T2;
    else if( theta_diff >= 0.5d0*PI .and. theta_diff < PI ) then
        angle = theta_r0 + T1 + T2;
    else if( theta_diff >= PI .and. theta_diff < 1.5d0*PI) then
        angle = theta_r0 - T1 - T2;
    else if( theta_diff >= 1.5d0*PI .and. theta_diff < 2d0*PI) then
        angle = theta_r0 - T1 + T2;
    end if

    if( angle < 0 ) angle = angle + 2d0 * PI;
    if( angle > 2d0 * PI ) angle = angle - 2d0 * PI;

    return;
end function asymptotic_angle
