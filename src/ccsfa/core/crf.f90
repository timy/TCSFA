
#include '../include/inc_plot_crf_track.h'
!#define CRF_NORM_ADD

    
    subroutine crf ( zs, hs, hm, dm, func, ds, ze, he, de, n )

        implicit none

        double complex:: a
        double complex:: cw
        double precision:: de
        double precision:: dm
        double precision:: ds
        double complex:: func
        double precision:: h
        double precision:: he
        double precision:: hm
        double precision:: hs
        integer:: i
        integer:: k
        integer:: n
        integer:: nr
        double complex:: u(7)
        double complex:: v
        double precision:: w(3)
        double precision:: w0
        double complex:: z(3)
        double complex:: z0
        double complex:: ze
        double complex:: zs


        u(1) = (  1.0d+00, 0.0d+00 )
        u(2) = (  0.8660254d+00, 0.5000000d+00 )
        u(3) = (  0.0000000d+00, 1.0000000d+00 )
        u(4) = (  0.9659258d+00, 0.2588190d+00 )
        u(5) = (  0.7071068d+00, 0.7071068d+00 )
        u(6) = (  0.2588190d+00, 0.9659258d+00 )
        u(7) = ( -0.2588190d+00, 0.9659258d+00 )

        h = hs;
        z0 = zs;
        n = 0;

        ! 1st calculation of ds.
        cw = func( z0 );


#ifdef CRF_PLOT_STEP
        open(CRF_PLOT_FILE_ID, file=CRF_PLOT_FILE_NAME);
#endif


#ifdef CRF_NORM_ADD
        w0 = dabs( dreal( cw ) ) + dabs ( dimag ( cw ) );
#else
        w0 = cdabs( cw );
#endif

        ds = w0; 

        ! if initial derivative reaches demanding precision
        if( w0 - dm <= 0 ) then
            ze = z0;
            he = h;
            de = w0;
            return;
        end if

        ! otherwise, start the downhill searching
        k = 1;
        i = 0;

        ! equilateral triangle walk pattern.
        v = ( -1d0, 0d0 );
        a = ( -0.5d0, 0.866d0 );

        do while ( w0 - dm > 0 )

            ! calculation of deviations w in the new test points.


            ! if necessary, record the searching path
#ifdef CRF_PLOT_STEP
            write(CRF_PLOT_FILE_ID, '(3(e15.8,2x))'), z0, w0;
#endif

            ! test point at the lower right
            z(1) = z0 + h * v * a;                    
            cw = func ( z(1) );

#ifdef CRF_NORM_ADD
            w(1) = dabs ( dreal( cw ) ) + dabs ( dimag ( cw ) );
#else
            w(1) = cdabs( cw );
#endif

            ! test point at the left
            z(2) = z0 + h * v;                        
            cw = func ( z(2) );

#ifdef CRF_NOMR_ADD
            w(2) = dabs ( dreal( cw ) ) + dabs ( dimag ( cw ) );
#else
            w(2) = cdabs( cw );
#endif

            ! test point at the upper right
            z(3) = z0 + h * conjg ( a ) * v;          
            cw = func ( z(3) );

#ifdef CRF_NORM_ADD
            w(3) = dabs ( dreal( cw ) ) + dabs ( dimag ( cw ) );
#else
            w(3) = cdabs( cw )
#endif

            n = n + 1;

            ! determination of w(nr) - the smallest of w(i).
            if( w(1) <= w(3) ) then
                if( w(1) < w(2) ) then
                    nr = 1;
                else
                    nr = 2;
                end if
            else
                if( w(2) < w(3) ) then
                    nr = 2;
                else
                    nr = 3;
                end if
            end if

            if( w0 - w(nr) >= 0 ) then 
                ! if the origin is not the lowest one compared with the test points
                ! set the lowest test point as the origin and continue the previous
                ! steps; the pushing direction and length of next step is altered
                k = 1;
                i = 0;
                a = ( 0.707d0, 0.707d0 );
                v = ( z(nr) - z0 ) / h;
                w0 = w(nr);
                z0 = z(nr);
                if( w0 - dm <= 0 ) then
                    exit;
                else
                    cycle;
                end if
            else
                ! if the origin is already the lowest one 
                ! compared with the test points
                if( k == 1 ) then
                    k = 2;
                    ! reduce the step length
                    if( h < hm ) then
                        ! fail to find a satisfying solution
                        exit;
                    else
                        ! reduce the step length, but don't 
                        ! change the direction
                        h = h * 0.25d0;
                        a = ( -0.5d0 , 0.866d0 );
                        cycle;
                    end if
                else if( k == 2) then
                    ! if after the reduction of length, 
                    ! there is still no solution
                    k = 3;
                    ! restore the step length
                    h = h * 4d0;
                    v = ( -1d0, 0d0 );
                    a = ( -0.5d0 , 0.866d0 );
                    cycle;
                else if( k == 3 ) then
                    i = i + 1;
                    ! rotate the walk pattern
                    if( i - 7 <= 0) then
                        ! if it hasn't tried enough times, keep 
                        ! trying changing the walk direction
                        v = u(i);
                        a = ( -0.5d0 , 0.866d0 );
                        cycle;
                    else
                        ! if merely changing direction doesn't work, 
                        ! reduce the step length and retry
                        if( h > hm ) then
                            h = h * 0.25d0;
                            i = 0;
                            v = ( -1d0, 0d0 );
                            a = ( -0.5d0 , 0.866d0 );
                            cycle
                        else
                            ! desperated >< !!
                            exit
                        end if
                    end if
                end if
            end if

        end do

        ze = z0;
        he = h;
        de = w0;

#ifdef CRF_PLOT_STEP
        close(CRF_PLOT_FILE_ID);
#endif


        return
    end subroutine crf
