#include '../include/inc_rk4.h'
#include '../include/inc_atom.h'
#include '../include/inc_misc.h'

subroutine rk4_re( t0_, tp_, x0_, vx0_, z0_, vz0_, f, ierr, W, px, pz, x, z, L, n_pass_x, n_pass_z )

    implicit none;
    integer, parameter:: nt = RK4_NT
    integer, parameter:: ne = RK4_NE
    integer, parameter:: nmax = RK4_NMAX
    double precision, intent(in):: t0_, tp_, x0_, vx0_, z0_, vz0_
    double precision:: eps = RK4_EPS
    external:: f;
    double precision:: t(nt), y(ne,nt), near_core_t0, near_core_t1
    double precision:: K1(ne), K2(ne), K3(ne), K4(ne);
    double precision:: h, p
    double precision:: ft, ftt, fy_old(ne), fy(ne, nmax), fyy(ne), fdyy(ne)
    integer:: it, jt, n_sub, i_next
    integer:: n_iter;

    integer, intent(out):: ierr
    double complex, intent(out):: W
    double precision, intent(out):: pz, px

    double complex:: w_near_core, dt
    double precision:: x, vx, z, vz, v2, r, energy, w_tail, L
    double precision:: Ex, Ez, ax, az, t_sub
    double precision:: charge, Ip
    integer:: i;
    double complex, external:: pulse_E_x, pulse_E_z

    integer:: n_pass_x, n_pass_z



#if MISC_PRINT > 3
    print*, 'rk4_re(): Initialization of RK4'
    print*, 't0_', t0_, 'tp_', tp_;
    print*, 'x0_', x0_, 'vx0_', vx0_;
    print*, 'z0_', z0_, 'vz0_', vz0_;
#endif

    h = ( tp_ - t0_ ) / ( nt - 1d0 );
    forall( i = 1:nt ) t(i) = t0_ + ( i - 1d0 ) * h;
    y(1,1) = x0_;
    y(2,1) = vx0_;
    y(3,1) = z0_
    y(4,1) = vz0_;



    W = (0d0, 0d0);
    n_pass_x = 0;
    n_pass_z = 0;

#ifdef MISC_PLOT
    call rk4_plot_init()
    call rk4_plot_open_file()
#endif

!    do it = 2, nt
    it = 1
    do while ( it .le. nt )

100     it = it + 1

        y(:,it) = 0d0;
        fy_old(:) = 0d0;

        n_sub = 2;
        h = ( t(it) - t(it-1) ) / ( n_sub - 1d0 );
        n_iter = 0;

        ! if the precision is not satisified, halve the step length and iterate
        do while( 1 )
            if ( n_sub > nmax ) then
                ! ierr = 2; ! too many steps still can't reach desired precision
                ! return;        

                t(it) = t(it-1) + 0.5 * ( t(it) - t(it-1) )
                it = it - 1
                go to 100 

            end if

            fy(:,1) = y(:,it-1)
            ft = t(it-1);

            do jt = 2, n_sub
                ftt = ft; 
                fyy(:) = fy(:,jt-1);
                call f( ne, ftt, fyy, fdyy );
                k1(:) = fdyy(:);

                ftt = ft + h / 2d0;
                fyy(:) = fy(:,jt-1) + h * k1(:) / 2d0;
                call f( ne, ftt, fyy, fdyy );
                k2(:) = fdyy(:);

                fyy(:) = fy(:,jt-1) + h * k2(:) / 2d0;
                call f( ne, ftt, fyy, fdyy );
                k3(:) = fdyy(:);

                ftt = ft + h;
                fyy(:) = fy(:,jt-1) + h * K3(:);
                call f( ne, ftt, fyy, fdyy );
                K4(:) = fdyy(:);

                fy(:,jt) = fy(:,jt-1) + &
                      h * (K1(:) + 2d0 * ( K2(:) + K3(:) ) + K4(:) ) / 6d0;

                ft = ft + h;
            end do

            ! at least two times of iterations are needed
            ! to verify the accuracy of the calculation

            p = 1d0;
            
            if( n_iter > 0 ) then

                p = maxval( dabs( fy_old(:) - fy(:, n_sub)) );

            end if

            fy_old(:) = fy(:,n_sub);

            if( p < eps ) then

#ifdef MISC_PLOT
                call rk4_plot_write( t(it-1), fy(1,n_sub), fy(2,n_sub), fy(3,n_sub), fy(4,n_sub) )
#endif
                
                charge = ATOM_CHARGE_Z;
                Ip = IONIZATION_IP;
                

                ! calculate the action W
                do i = 2, n_sub
                    
                    if( fy(1,i-1) * fy(1,i) < 0 ) n_pass_x = n_pass_x + 1;
                    if( fy(3,i-1) * fy(3,i) < 0 ) n_pass_z = n_pass_z + 1;

                    x = ( fy(1,i-1) + fy(1,i) ) / 2d0;
                    vx = ( fy(2,i-1) + fy(2,i) ) / 2d0;
                    z = ( fy(3,i-1) + fy(3,i) ) / 2d0;
                    vz = ( fy(4,i-1) + fy(4,i) ) / 2d0;
                    v2 = vz*vz + vx*vx;
                    r = dsqrt( x*x + z*z );

                    ! if the electron gets too close to the core
                    if( r < RK4_R_THRESHOLD ) then
                        
                        !ierr = 4; ! high-order term we cannot handle currently
                        !return;
                        
                        print*, 'Rk4 threshold activated!'

                        call kepler_near_core( z, x, vz, vx, charge, Ip, dt, &
                              y(3,it), y(1,it), y(4,it), y(2,it), w_near_core )
                        w = w + w_near_core
                        ! get current time
                        near_core_t0 = t(it-1) + ( i - 1.5d0 ) * h
                        ! get the time at the exit
                        near_core_t1 = near_core_t0 + dreal( dt )
                        ! get the index of the time at the exit
                        i_next = ceiling( near_core_t1 / ( ( tp_ - t0_ ) / ( nt - 1d0 ) ) )
                        if( i_next > nt ) then ! currently we don't consider the near-core for pulse-off ..
                            ierr = 4
                            return
                        end if
                        ! reset the element of the time array
                        t(i_next) = near_core_t1
                        ! reset the loop index ( skip the intermediate time steps )
                        it = i_next
                        
                        exit

                    end if

                    if( REPRESENTATION_OPT == 'S' ) then

                        energy = 0.5d0 * v2 - charge / r;
                        ! here modified for urgent use of screened potential !!!
                        !energy = 0.5d0 * v2 - 1d0 / r - (charge - 1d0) * dexp(-5d0*r) / r

                    elseif( REPRESENTATION_OPT == 'W' ) then

                        t_sub = t(it-1) + (i-1.5d0) * h;
                        Ex = dreal( pulse_E_x( dcmplx(t_sub) ) );
                        Ez = dreal( pulse_E_z( dcmplx(t_sub) ) );
                        ax = -Ex - charge * x / r**3;
                        az = -Ez - charge * z / r**3;
                        energy = 0.5d0 * v2 - charge / r + Ex * x + Ez * z + ax * x + az * z;

                    end if

                    W = W + dcmplx( ( energy + Ip ) * h, 0d0 );
                end do

                exit

            else

                h = h / 2d0;
                n_sub = n_sub + n_sub - 1;
                n_iter = n_iter + 1;

            end if

        end do
        ! we have finished one desired interval in the array,
        ! then go to the next interval

#if MISC_PRINT > 4
        write(*,'(a,i6,a,i6)'), 'number of iteration for sub-step ', &
              it, ': ', n_iter - 1 ;
#endif
        y(:,it) = fy_old(:);

    end do

    x = y(1,nt);
    vx = y(2,nt);
    z = y(3,nt);
    vz = y(4,nt);

    ! calculate asymptotic momentum

    call coulomb_tail( z, x, vz, vx, charge, L, pz, px, w_tail, ierr )
    ! energy is less than 0: ierr == 3
    if( ierr > 0 ) then 
        return
    end if

    ! W_cen for the W-representation
    if( REPRESENTATION_OPT == 'W' ) then
        w = w + dcmplx(w_tail, 0d0);
    end if

#ifdef MISC_PLOT
    call rk4_plot_close_file()
#endif
 
    ierr = 0

    return;
end subroutine rk4_re
