#define PULSE_NT 2000
#define PULSE_FILE_ID 201
#define PULSE_FILE_NAME 'dat/pulse.dat'


module mod_pulse

    implicit none;

    double precision:: E0;
    double precision:: om;
    double precision:: nc;
    double precision:: xi;
    double precision:: ph;
    double precision, parameter::  PI = 2d0*asin(1d0);
    
    double precision, parameter:: n_head = 0.01d0
    double precision, parameter:: n_tail = 0.01d0
    double precision, parameter:: n_const = 8.98d0
    double precision:: head
    double precision:: tail
    double precision:: period
    double precision:: tp;

    ! t0 is used only when calculating alpha, actually it has no influence whatever value
    ! it has, since it vanishes in the calculation of x0 and z0
    double complex:: t0 = ( 0d0, 0d0 );     

contains

    subroutine set_pulse( E0_, om_, nc_, xi_, ph_ )

        implicit none;
        double precision, intent(in):: E0_, om_, nc_, xi_, ph_

        E0 = E0_;
        om = om_;
        nc = nc_;
        xi = xi_;
        ph = ph_;
        ! here should be modified for different kinds of problems
        ! here is for ramping
        ! period = 2d0 * PI / om
!         head = period * n_head
!         tail = period * (n_head + n_const)
!         tp = period * ( n_head + n_const + n_tail )
        ! ramping end

        ! for general sin2 or const
        tp = 2d0 * PI * nc / om;

        return;
    end subroutine set_pulse


    subroutine set_pulse_t0( t0_ )

        implicit none;
        double complex, intent(in):: t0_

        t0 = t0_;

        return;
    end subroutine set_pulse_t0

end module mod_pulse



! double complex function pulse_A_z(t) result(A_z)
!     implicit none
!     double complex, intent(in):: t
!     double complex, external:: pulse_A_z_sin2
!     A_z = pulse_A_z_sin2(t)
!     return
! end function pulse_A_z
    
! double complex function pulse_A_x(t) result(A_x)
!     implicit none
!     double complex, intent(in):: t
!     double complex, external:: pulse_A_x_sin2
!     A_x = pulse_A_x_sin2(t)
!     return
! end function pulse_A_x

! double complex function pulse_E_z(t) result(E_z)
!     implicit none
!     double complex, intent(in):: t
!     double complex, external:: pulse_E_z_sin2
!     E_z = pulse_E_z_sin2(t)
!     return
! end function pulse_E_z

! double complex function pulse_E_x(t) result(E_x)
!     implicit none
!     double complex, intent(in):: t
!     double complex, external:: pulse_E_x_sin2
!     E_x = pulse_E_x_sin2(t)
!     return
! end function pulse_E_x

! double complex function pulse_alpha_z(t) result(alpha_z)
!     implicit none
!     double complex, intent(in):: t
!     double complex, external:: pulse_alpha_z_sin2
!     alpha_z = pulse_alpha_z_sin2(t)
!     return
! end function pulse_alpha_z

! double complex function pulse_alpha_x(t) result(alpha_x)
!     implicit none
!     double complex, intent(in):: t
!     double complex, external:: pulse_alpha_x_sin2
!     alpha_x = pulse_alpha_x_sin2(t)
!     return
! end function pulse_alpha_x

! ! double complex function pulse_alpha2_z(t) result(alpha2_z)
! !     implicit none
! !     double complex, intent(in):: t
! !     double complex, external:: pulse_alpha2_z_sin2
! !     alpha2_z = pulse_alpha2_z_sin2(t)
! !     return
! ! end function pulse_alpha2_z

! ! double complex function pulse_alpha2_x(t) result(alpha2_x)
! !     implicit none
! !     double complex, intent(in):: t
! !     double complex, external:: pulse_alpha2_x_sin2
! !     alpha2_x = pulse_alpha2_x_sin2(t)
! !     return
! ! end function pulse_alpha2_x


double precision function get_tp()
    use mod_pulse, only: tp;
    implicit none;

    get_tp = tp;

    return;
end function get_tp




subroutine plot_pulse()
    use mod_pulse, only: tp
    implicit none
    double complex, external:: pulse_A_z, pulse_A_x, pulse_E_z, pulse_E_x, &
          pulse_alpha_z, pulse_alpha_x, pulse_alpha2_z, pulse_alpha2_x
    double precision:: dt
    integer:: i
    double complex:: t

    open( PULSE_FILE_ID, file=PULSE_FILE_NAME );
    dt = tp / ( PULSE_NT - 1d0 );
    
    do i = 1, PULSE_NT

        t = dcmplx( (i - 1d0) * dt, 0d0 );
        write( PULSE_FILE_ID, '(7(e15.8,1x))'), dreal(t), &
              dreal(pulse_A_z(t)), dreal(pulse_A_x(t)), &
              dreal(pulse_E_z(t)), dreal(pulse_E_x(t)), &
              dreal(pulse_alpha_z(t)), dreal(pulse_alpha_x(t))!, &
         !     dreal(pulse_alpha2_z(t)), dreal(pulse_alpha2_x(t));
        
    end do

    close( PULSE_FILE_ID );

end subroutine plot_pulse
