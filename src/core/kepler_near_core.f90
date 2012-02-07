subroutine kepler_near_core( x0, y0, vx0, vy0, Z, Ip, dt, x1, y1, vx1, vy1, w )
    
    implicit none
    double precision, intent(in):: x0, y0, vx0, vy0, Z, Ip
    double precision, intent(out):: x1, y1, vx1, vy1
    double complex, intent(out):: w, dt
    double precision:: r0, energy, epsilon, L, k, phi0
    double precision, parameter:: pi = 2d0 * asin(1d0)
    double complex, parameter:: eye = dcmplx( 0d0, 1d0 )
    double complex:: a1, a2, theta0, gamma, vx_polar, vy_polar, theta_v_income_polar, theta_v_outcome_polar, &
          theta_v_outcome
    double precision:: a, b

    r0 = dsqrt( x0 * x0 + y0 * y0 )

    energy = 0.5 * ( vx0 * vx0 + vy0 * vy0 ) - Z / r0
    
    L = x0 * vy0 - y0 * vx0

    k = L * L / Z

    epsilon = dsqrt( 1d0 + 2d0 * energy * L * L / ( Z * Z ) )

    phi0 = datan2( y0, x0 )

    theta0 = cdacos( dcmplx( ( k - r0 ) / ( r0 * epsilon ) ) )
    
    if( ( x0 * vx0 + y0 * vy0 ) * L .lt. 0d0 ) theta0 = - theta0

    gamma = phi0 - theta0
    
    !print*, 180 * phi0 / pi, 180 * theta0 / pi, 180 * gamma / pi


    ! calculate the time difference

    a1 = cdatanh( ( epsilon - 1d0 ) / cdsqrt( dcmplx( epsilon*epsilon - 1d0 ) ) * cdtan( 0.5d0 * theta0 ) )
    a1 = 2d0 * a1 / dcmplx( epsilon * epsilon - 1d0 )**1.5
    a2 = - epsilon * cdsin( theta0 ) / ( ( epsilon * epsilon - 1d0 ) * ( 1d0 + epsilon * cdcos( theta0 ) ) )
    dt = 2d0 * k * k / L * ( a1 + a2 )
    
    ! calculate the new position x1, y1
    
    x1 = dreal( r0 * cdcos( gamma - theta0 ) )
    y1 = dreal( r0 * cdsin( gamma - theta0 ) )

    ! calculate the new velocity vx1, vy1
    
    call inv_rot_matrix( dcmplx(vx0), dcmplx(vy0), vx_polar, vy_polar )
    theta_v_income_polar = cdacos( vx_polar / cdsqrt( vx_polar * vx_polar + vy_polar * vy_polar ) )
    if( dreal( vy_polar ) < 0d0 ) theta_v_income_polar = - theta_v_income_polar
    theta_v_outcome_polar = pi - theta_v_income_polar
    theta_v_outcome = theta_v_outcome_polar + gamma
    
    !print*, 180 * theta_v_income_polar / pi, 180 * theta_v_outcome_polar / pi, 180 * theta_v_outcome / pi

    vx1 = dreal( dsqrt( vx0 * vx0 + vy0 * vy0 ) * cdcos( theta_v_outcome ) )
    vy1 = dreal( dsqrt( vx0 * vx0 + vy0 * vy0 ) * cdsin( theta_v_outcome ) )
   
    ! calculate the action W

    a = 0.5 * L * L / ( k * k )
    b = Z / k
    a1 = cdatanh( ( epsilon - 1d0 ) * cdtan( 0.5d0 * theta0 ) / cdsqrt( dcmplx( epsilon * epsilon - 1d0 ) ) )
    a1 = 2d0 * ( a + Ip - a * epsilon * epsilon + b * ( epsilon * epsilon - 1d0 ) ) * a1
    a1 = a1 / dcmplx( epsilon * epsilon - 1d0 )**1.5
    a2 = -epsilon * ( Ip + a * ( epsilon * epsilon - 1d0 ) ) * cdsin( theta0 )
    a2 = a2 / ( ( epsilon * epsilon - 1d0 ) * ( 1d0 + epsilon * cdcos( theta0 ) ) )
    w = 2d0 * k * k * ( a1 + a2 ) / L
    
    return

! ============================================================
contains

    ! ----------------------------------------
    subroutine rot_matrix( x_old, y_old, x_new, y_new )

        implicit none

        double complex:: x_old, y_old, x_new, y_new

        x_new = cdcos( gamma ) * x_old - cdsin( gamma ) * y_old
        y_new = cdsin( gamma ) * x_old + cdcos( gamma ) * y_old
        
        return
    end subroutine rot_matrix

    ! ----------------------------------------
    subroutine inv_rot_matrix( x_old, y_old, x_new, y_new )

        implicit none

        double complex:: x_old, y_old, x_new, y_new
        
        x_new = cdcos( gamma ) * x_old + cdsin( gamma ) * y_old
        y_new = -cdsin( gamma ) * x_old + cdcos( gamma ) * y_old
        
        return
        
    end subroutine inv_rot_matrix

    ! ----------------------------------------
    double complex function cdatanh( z )

        implicit none

        double complex:: z
        
        cdatanh = 0.5d0 * ( cdlog( 1d0 + z ) - cdlog( 1d0 - z ) )
        
        return
        
    end function cdatanh
    
    ! ----------------------------------------
    double complex function cdacos( z )
        
        implicit none
        
        double complex:: z
        
        cdacos = 0.5 * pi + eye * cdlog( eye * z + cdsqrt( 1d0 - z * z ) )
        
    end function cdacos

end subroutine kepler_near_core
