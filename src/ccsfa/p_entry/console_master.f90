#include "../include/inc_grid_size.h"
#include "../include/inc_grid_boundary.h"

subroutine console_master( mpi_info )

    use mod_mmff_mpi_info, only: T_mmff_mpi_info;

    implicit none;
    type(T_mmff_mpi_info), intent(in):: mpi_info;
    double precision, pointer:: px_0_(:), pz_0_(:), x_0_(:), z_0_(:), px_inf_(:), pz_inf_(:)
    double complex, pointer:: amp_M_(:), ts_(:)
    integer, pointer:: n_ts_in_p0_(:), n_pass_x_(:), n_pass_z_(:), ierr_(:)
    integer, pointer:: a_task(:), n_ts(:)
    integer, allocatable:: n_p0_in_proc(:), n_ts_in_proc(:)
    external:: null_dbl, null_int, null_dcp
    integer:: i_pz, i_px, first, last, i_proc, n_proc, n_all_ts, count, i_ts
    integer:: h_q1, h_q2, h_q3
    integer:: h_a1, h_a2, h_a3, h_a4, h_a5, h_a6, h_a7, h_a8, h_a9, h_a10, h_a11, h_a12, h_a13, h_a14
    integer:: n_data
    integer:: time_start(6)

    interface
        include 'mmff_console_interface.f90'
    end interface

    ! initialize
    call mmff_init( mpi_info );

    !* send px_0 and pz_0
    n_data = N_P0;

    call mmff_reset( n_data, 0, h_q1 );
    call mmff_create_data_dbl( px_0_, null_dbl, 'px_0', h_a1 );
    call mmff_create_data_dbl( pz_0_, null_dbl, 'pz_0', h_a2 );
    !! initialize the parameters
    count = 0;
    do i_pz = 1, N_PZ
        do i_px = 1, N_PX
            count = count + 1;
            px_0_(count) = PX_BEGIN + ( i_px - 1d0 ) * D_PX;
            pz_0_(count) = PZ_BEGIN + ( i_pz - 1d0 ) * D_PZ;
        end do
    end do
    !!
    call mmff_send_data( h_q1 );
    call mmff_delete_data_dbl( px_0_, h_a1 );
    call mmff_delete_data_dbl( pz_0_, h_a2 );


    n_proc = mpi_info % n_proc - 1;
    allocate( n_p0_in_proc( n_proc ) );

    !* recv n_ts
    call  timer_tic( time_start )
    n_data = N_P0;

    call mmff_reset( n_data, 0, h_q2 );
    !! get how many p0 for each proc
    call mmff_get_task( h_q2, a_task );
    forall(i_proc=1:n_proc) n_p0_in_proc(i_proc) = a_task(i_proc); 
    call mmff_set_task( h_q2, a_task );
    !!
    call mmff_create_data_int( n_ts_in_p0_, null_int, 'n_ts', h_a3 );
    call mmff_recv_data( h_q2 );
    
    call  timer_toc_with_log( time_start )


    !* get how many traj for each proc
    allocate( n_ts_in_proc(n_proc) );
    do i_proc = 1, n_proc !loop over all processes
        first = sum( n_p0_in_proc(1:i_proc-1) ) + 1; !
        last = sum( n_p0_in_proc(1:i_proc) );
        n_ts_in_proc(i_proc) = sum( n_ts_in_p0_( first:last ) );
    end do
    n_all_ts = sum( n_ts_in_proc );



    !* broadcast the task of each proc to slaves
    call mmff_broadcast_master_int( n_proc, n_ts_in_proc );


    !* recv W, p0_inf
    n_data = n_all_ts; 
    call mmff_reset( n_data, 0, h_q3 );
    call mmff_get_task( h_q3, a_task );
    forall(i_proc=1:n_proc) a_task(i_proc) = &
          n_ts_in_proc(i_proc);
    call mmff_set_task( h_q3, a_task );
    call mmff_create_data_dbl( px_0_, null_dbl, 'px_0', h_a4 );
    call mmff_create_data_dbl( pz_0_, null_dbl, 'pz_0', h_a5 );
    call mmff_create_data_dcp( ts_, null_dcp, 'ts', h_a6 );
    call mmff_create_data_dbl( x_0_, null_dbl, 'x_t0', h_a7 );
    call mmff_create_data_dbl( z_0_, null_dbl, 'z_t0', h_a8);
    call mmff_create_data_dbl( px_inf_, null_dbl, 'px_inf', h_a9 );
    call mmff_create_data_dbl( pz_inf_, null_dbl, 'pz_inf', h_a10 );
    call mmff_create_data_dcp( amp_M_, null_dcp, 'amp_M', h_a11 );
    call mmff_create_data_int( n_pass_x_, null_int, 'n_pass_x', h_a12);
    call mmff_create_data_int( n_pass_z_, null_int, 'n_pass_z', h_a13);
    call mmff_create_data_int( ierr_, null_int, 'ierr', h_a14 );
    call mmff_recv_data( h_q3 );

    !! data analysis and output
    open( 222, file = 'dat/data.dat' );
    do i_ts = 1, n_all_ts
        write( 222, '(10(e15.8, 1x), 3(i2,1x))' ), &
              px_0_(i_ts), pz_0_(i_ts), &
              ts_(i_ts), &
              x_0_(i_ts), z_0_(i_ts), &
              px_inf_(i_ts), pz_inf_(i_ts), &
              amp_M_(i_ts), &
              n_pass_x_(i_ts), n_pass_z_(i_ts), &
              ierr_(i_ts);
    end do
    close( 222 );

    !* finalization
    deallocate( n_p0_in_proc );
    deallocate( n_ts_in_proc );
!!$
    call mmff_delete_data_int( n_ts_in_p0_, h_a3 );
    call mmff_delete_data_dbl( px_0_, h_a4 );
    call mmff_delete_data_dbl( pz_0_, h_a5 );
    call mmff_delete_data_dcp( ts_, h_a6 );
    call mmff_delete_data_dbl( x_0_, h_a7 );
    call mmff_delete_data_dbl( z_0_, h_a8)
    call mmff_delete_data_dbl( px_inf_, h_a9 );
    call mmff_delete_data_dbl( pz_inf_, h_a10 );
    call mmff_delete_data_dcp( amp_M_, h_a11 );
    call mmff_delete_data_int( n_pass_x_, h_a12 );
    call mmff_delete_data_int( n_pass_z_, h_a13 );
    call mmff_delete_data_int( ierr_, h_a14 );

    call mmff_final();

    print*, 'Done!', n_all_ts, 'traj were output'
    return;

end subroutine console_master
