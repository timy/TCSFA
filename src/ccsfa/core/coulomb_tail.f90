subroutine coulomb_tail( x, y, vx, vy, charge, L, px, py, w, ierr )
    
    implicit none;
    double precision, intent(in):: x, y, vx, vy, charge
    double precision, intent(out):: L, px, py, w
    integer, intent(out):: ierr
    double precision:: r0, E, asym_ang, eps, phi0, k, a, b, c, th_0, th_inf
    double precision, external:: angle, to_deg, renorm_ang
    double precision, parameter:: PI = 2d0 * asin( 1d0 )

    ! if there is no long range central force
    if( dabs(charge) < 1d-8 ) then
        px = vx;
        py = vy;
        w = 0d0;
        return;
    end if

    r0 = dsqrt( x*x + y*y );
    E = - charge / r0 + 0.5d0 * ( vx*vx + vy*vy );
    if( E <= 0 ) then
        ierr = 2; 
        return;
    end if

    L = x*vy - y*vx;
    eps = dsqrt( 1d0 + 2d0 * E * L*L / ( charge*charge ) );
    phi0 = datan2( y, x );
!!$    k = L * L / charge;
!!$    a = k / ( eps*eps - 1d0 );
!!$    b = a * dsqrt( eps*eps - 1d0 );
!!$    c = eps * a;

    th_inf = dacos( -1d0 / eps );
    if( L < 0 ) th_inf = -th_inf;
    th_0 = dacos( ( L*L / (charge*r0) - 1d0 ) / eps );
    if( ( L < 0 .and. angle( x, y, vx, vy ) <= PI/2d0) .or. &
          ( L >= 0 .and. angle( x, y, vx, vy ) > PI/2d0)) then
        th_0 = -th_0;
    end if
    asym_ang =  phi0 + th_inf - th_0;
    px = dsqrt( 2d0 * E ) * dcos( asym_ang );
    py = dsqrt( 2d0 * E ) * dsin( asym_ang );

    w = ( 2d0 * L / dsqrt( eps*eps - 1d0 ) ) * &
          datanh( dsqrt( ( eps - 1d0 ) / ( eps + 1d0 ) ) * &
          dtan( th_0 / 2d0 ) );
    ierr = 0;

!!$    write(*,'(6(a,f15.8,2x))'), 'E=',E, 'L=',L, 'eps=',eps, 'px=',px, 'py=',py;
!!$    write(*,'(3(a,f15.8,2x))'), 'phi0=', to_deg(phi0), 'th_0=', to_deg(th_0), &
!!$          'th_inf=', to_deg(th_inf), 'alpha=', to_deg(asym_ang);

    return;
    
end subroutine coulomb_tail

double precision function angle( x1, y1, x2, y2 )
    implicit none
    double precision, intent(in):: x1, y1, x2, y2
    angle = dacos( ( x1*x2 + y1*y2 ) / dsqrt((x1*x1+y1*y1)*(x2*x2+y2*y2)) );
    return;
end function angle

double precision function to_deg( ang )
    implicit none
    double precision, parameter:: PI = 2d0*dasin(1d0)
    double precision, intent(in):: ang
    to_deg = ang * 180 / PI;
    return;
end function to_deg

double precision function renorm_ang( ang )
    implicit none
    double precision, parameter:: PI = 2d0*dasin(1d0)
    double precision, intent(in):: ang
    double precision:: n_cyc
    n_cyc = ang / (2d0*PI);
    renorm_ang = ang - 2d0*PI*floor(n_cyc);
    return;
end function renorm_ang


!!$program main
!!$
!!$    implicit none
!!$
!!$    double precision:: z, x, vz, vx, charge, asym_ang, pz, px, w
!!$    integer:: ierr
!!$   
!!$    z = 2.03d0;
!!$    x = -2.70d0;
!!$    vz = -0.8259d0;
!!$    vx = -0.6312d0;
!!$    charge = 1d0;
!!$    
!!$    call coulomb_tail( z, x, vz, vx, charge, asym_ang, pz, px, w, ierr );
!!$
!!$    write(*, '(3(a,f15.8,2x))'), 'pz=', pz, 'px=', px, 'w', w;
!!$
!!$end program main
