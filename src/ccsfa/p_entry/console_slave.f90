#include "../include/inc_grid_size.h"
#include "../include/inc_field.h"
#include "../include/inc_atom.h"
#include "../include/inc_ts_guess.h"

subroutine console_slave( mpi_info )

    use mod_mmff_mpi_info, only: T_mmff_mpi_info;
    use mod_pulse, only: set_pulse, set_pulse_t0

    implicit none;
    ! MPI_related
    type(T_mmff_mpi_info), intent(in):: mpi_info;
    double precision, pointer:: px_0_(:), pz_0_(:), x_0_(:), z_0_(:), px_inf_(:), pz_inf_(:)
    double complex, pointer:: amp_M_(:), ts_(:)
    integer, pointer:: n_ts_in_p0_(:), n_ts_in_proc(:), ierr_(:), n_pass_x_(:), n_pass_z_(:)
    integer, pointer:: a_task(:), n_ts(:)
    double precision, allocatable:: px_0(:), pz_0(:)
    external:: null_dbl, null_int, null_dcp
    integer:: i_proc, i_p0, i_ts, count
    integer:: h_q1, h_q2, h_q3
    integer:: h_a1, h_a2, h_a3, h_a4, h_a5, h_a6, h_a7, h_a8, h_a9, h_a10, h_a11, h_a12, h_a13, h_a14
    integer:: n_data, n_proc, n_p0, n_all_ts

    ! CCSFA_related
    double complex, allocatable:: a_ts_guess(:,:)
    double precision, external:: get_tp 
    double complex, external:: SPE
    double precision:: tp



    interface
        include 'mmff_console_interface.f90'
    end interface


    ! initialize
    call mmff_init( mpi_info );

    !* recv px_0 and pz_0
    n_data = N_P0;
    call mmff_reset( n_data, 0, h_q1 );
    call mmff_create_data_dbl( px_0_, null_dbl, 'px_0', h_a1 );
    call mmff_create_data_dbl( pz_0_, null_dbl, 'pz_0', h_a2 );
    call mmff_recv_data( h_q1 );
    n_p0 = n_data;
    allocate( px_0(n_p0) );
    allocate( pz_0(n_p0) );
    do i_p0 = 1, n_p0
        px_0(i_p0) = px_0_(i_p0);
        pz_0(i_p0) = pz_0_(i_p0);
    end do
    call mmff_delete_data_dbl( px_0_, h_a1 );
    call mmff_delete_data_dbl( pz_0_, h_a2 );


    n_proc = mpi_info % n_proc - 1;


    !* send n_ts
    n_data = N_P0;
    call mmff_reset( n_data, 0, h_q2 );
    call mmff_create_data_int( n_ts_in_p0_, null_int, 'n_ts', h_a3 );
    !! -call local_minima to estimate the number of ts
    ! ============================================================


    call set_pulse( E0, OM, NC, XI );
    call set_pulse_t0( (1d0, 0d0) );

    allocate( a_ts_guess(n_p0, LMS_MAX_COUNT) );
    do i_p0 = 1, n_p0

        call set_p0( px_0(i_p0), pz_0(i_p0) );
        tp = get_tp();
        call local_minima( LMS_RE_LOWER, LMS_RE_UPPER, &
              LMS_IM_LOWER, LMS_IM_UPPER, SPE, &
              n_ts_in_p0_(i_p0), a_ts_guess(i_p0,:) );

    end do
    ! ============================================================
    !!
    call mmff_send_data( h_q2 );


    !* broadcast the task of each proc to slaves
    call mmff_broadcast_slave_int( n_ts_in_proc );
    n_all_ts = sum( n_ts_in_proc );

    !* send W, p0_inf
    n_data = n_all_ts;
    call mmff_reset( n_data, 0, h_q3 );
    call mmff_get_task( h_q3, a_task );
    forall(i_proc=1:n_proc) a_task(i_proc) = &
          n_ts_in_proc(i_proc);
    call mmff_set_task( h_q3, a_task );
    call mmff_create_data_dbl( px_0_, null_dbl, 'px_0', h_a4 );
    call mmff_create_data_dbl( pz_0_, null_dbl, 'pz_0', h_a5 );
    call mmff_create_data_dcp( ts_, null_dcp, 'ts', h_a6 );
    call mmff_create_data_dbl( x_0_, null_dbl, 'x_tp', h_a7 );
    call mmff_create_data_dbl( z_0_, null_dbl, 'z_tp', h_a8 );
    call mmff_create_data_dbl( px_inf_, null_dbl, 'px_inf', h_a9 );
    call mmff_create_data_dbl( pz_inf_, null_dbl, 'pz_inf', h_a10 );
    call mmff_create_data_dcp( amp_M_, null_dcp, 'W', h_a11 );
    call mmff_create_data_int( n_pass_x_, null_int, 'n_pass_x', h_a12);
    call mmff_create_data_int( n_pass_z_, null_int, 'n_pass_z', h_a13);
    call mmff_create_data_int( ierr_, null_int, 'ierr', h_a14);

    !! calculate px_inf, pz_inf, action_w
    ! ============================================================
    count = 0;
    do i_p0 = 1, n_p0
        do i_ts = 1, n_ts_in_p0_(i_p0);
            count = count + 1;
            px_0_(count) = px_0(i_p0);
            pz_0_(count) = pz_0(i_p0);

            call propagate_with_single_p0( &
                  px_0(i_p0), pz_0(i_p0), &
                  a_ts_guess(i_p0, i_ts), &
                  ts_(count), amp_M_(count), &
                  x_0_(count), z_0_(count), &
                  px_inf_(count), pz_inf_(count), &
                  n_pass_x_(count), n_pass_z_(count), &
                  ierr_(count) );
!            print*, 'rank', mpi_info%i_rank, 'n_ts', n_ts_in_p0_(i_p0), &
!                  'px_0', px_0(i_p0), 'pz_0', pz_0(i_p0)
        end do
    end do
    ! ============================================================
    !!
    call mmff_send_data( h_q3 );
!!$    !* finalization
    deallocate( px_0 );
    deallocate( pz_0 );
    deallocate( n_ts_in_proc );
    deallocate( a_ts_guess );
!!$
!!$
    call mmff_delete_data_int( n_ts_in_p0_, h_a3 );
    
    call mmff_delete_data_dbl( px_0_, h_a4 );
    call mmff_delete_data_dbl( pz_0_, h_a5 );
    call mmff_delete_data_dcp( ts_, h_a6 );
    call mmff_delete_data_dbl( x_0_, h_a7 );
    call mmff_delete_data_dbl( z_0_, h_a8 );
    call mmff_delete_data_dbl( px_inf_, h_a9 );
    call mmff_delete_data_dbl( pz_inf_, h_a10 );
    call mmff_delete_data_dcp( amp_M_, h_a11 );
    call mmff_delete_data_int( n_pass_x_, h_a12 )
    call mmff_delete_data_int( n_pass_z_, h_a13 )
    call mmff_delete_data_int( ierr_, h_a14 );

    call mmff_final();


    return;

end subroutine console_slave


