!------------------------------------------------------------------
!    Double Exponential Numerical Integration for Complex Type 
!------------------------------------------------------------------

double complex function cedint( f, ni, a, b )
        implicit none;
        double precision, parameter:: pi4 = 0.78539816339744830962d0
        ! --------------------------------
        integer:: n, ni;
        double complex, intent(in):: a, b
        double complex, external:: f;
        ! ---------------------------------
        double complex:: h, ss, z, exz, hcos, hsin, s, dxdz, x, w;
        integer:: k;
        
        n = ni;
        n = n / 2;
        h = 5d0 / n;
        
        ss = 0d0;
        do k = -n, n
                z = h * k;
                exz = cdexp( z );
                hcos = exz + 1d0 / exz;
                hsin = exz - 1d0 / exz;
                s = cdexp( pi4 * hsin )
                w = s + 1d0 / s;
                x = ( b * s + a / s ) / w;
                if ( x /= a .and. x /= b) then
                        dxdz = hcos/ ( w * w );
                        ss = ss + f(x) * dxdz;
                end if
        end do

        do k = -n, n
                z = h * ( dcmplx(k) + 0.5d0 );
                exz = cdexp(z);
                hcos = exz + 1d0 / exz;
                hsin = exz - 1d0 / exz;
                s = cdexp( pi4 * hsin );
                w = s + 1d0 / s;
                x = ( b*s + a/s ) / w;
                if ( x /= a .and. x /= b ) then
                        dxdz = hcos / ( w * w );
                        ss = ss + f(x) * dxdz;
                end if
        end do
        cedint = h * ( b - a ) * pi4 * ss;
        
        return;
end function cedint
