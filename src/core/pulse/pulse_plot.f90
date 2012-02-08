#include '../../include/inc_field.h'

subroutine pulse_plot
    implicit none
    integer, parameter:: n_samp = 200
    double complex, external:: PULSE_A_Z, PULSE_A_X
    double complex, external:: PULSE_E_Z, PULSE_E_X
    double complex, external:: PULSE_ALPHA_Z, PULSE_ALPHA_X
    double complex, external:: PULSE_ALPHA2_Z, PULSE_ALPHA2_X
    double complex:: dt, t
    double precision:: A_z, A_x, E_z, E_x, alpha_z, alpha_x, alpha2_z, alpha2_x

    dt = (Tp -  0d0) / (n_samp - 1d0)
    open(FID_PULSE, file=FNM_PULSE)
    do i = 1, n_samp
        t = dcmplx( 0d0+(i - 1d0)*dt, 0d0 )
        A_z = dreal( PULSE_A_Z( t ) )
        A_x = dreal( PULSE_A_X( t ) )
        E_z = dreal( PULSE_E_Z( t ) )
        E_x = dreal( PULSE_E_X( t ) )
        alpha_z = dreal( PULSE_ALPHA_Z( t ) )
        alpha_x = dreal( PULSE_ALPHA_X( t ) )
        alpha2_z = dreal( PULSE_ALPHA2_Z( t ) )
        alpha2_x = dreal( PULSE_ALPHA2_X( t ) )
        write(FID_PULSE, '(8(f15.8,1x))'), A_z, A_x, E_z, E_x, &
             alpha_z, alpha_x, alpha2_z, alpha2_x
    end do
    close(FID_PULSE)
end subroutine pulse_plot
