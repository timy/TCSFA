#include '../include/inc_field.h'
#include '../include/inc_ts_guess.h'
#include '../include/inc_grid_boundary.h'
#include '../include/inc_grid_size.h'
#include '../include/inc_rk4.h'

subroutine slave_thread( n_node, rank )
    implicit none

    integer, intent(in):: n_node, rank
    character(len=64):: file_name

    integer, external:: amount_of_job, start_of_job
    integer:: n_p0, p0_start
    integer:: ix, iz, i_p0, i_ts, index
    double precision, allocatable:: px(:), pz(:)
    integer, allocatable:: n_ts(:)
    double complex, allocatable:: ts_guess(:,:)
    double precision:: rand_px_0, rand_pz_0

    double complex:: ts, Mp
    double precision:: x0, z0, px_inf, pz_inf, L, err_spe
    integer:: n_pass_x, n_pass_z, ierr, n_step
    integer:: seed(1)
    integer, parameter:: tag = 0

    double complex, external:: SPE

    ! organize my tasks
    n_p0 = amount_of_job( N_P0, n_node, rank )
    p0_start = start_of_job( N_P0, n_node, rank )

    ! output info
    write( file_name, '(a,i3,a)' ), 'dat/info_', rank, '.dat'
    open( 101, file=file_name )
    call write_time( 101 )
    

    ! init the momentum according to my rank
    allocate( px(n_p0) )
    allocate( pz(n_p0) )
    seed(1) = rank
    call random_seed( put=seed )
    do i_p0 = p0_start, p0_start+n_p0-1
        ix = i_p0 / N_PZ
        iz = mod( i_p0, N_PZ )
        index = i_p0 - p0_start + 1
        call random_number( harvest = rand_px_0 )
        call random_number( harvest = rand_pz_0 )
        !rand_px_0 = 0.5d0
        !rand_pz_0 = 0.5d0
!!$        px(index) = PX_BEGIN + ( ix + rand_px_0 - 0.5d0 ) * D_PX
!!$        pz(index) = PZ_BEGIN + ( iz + rand_pz_0 - 0.5d0 ) * D_PZ

        ! completely random distribution
        px(index) = PX_BEGIN + (PX_END - PX_BEGIN) * rand_px_0
        pz(index) = PZ_BEGIN + (PZ_END - PZ_BEGIN) * rand_pz_0

    end do

    ! output px, pz
    write( file_name, '(a,i3,a)' ), 'dat/p_', rank, '.dat'
    open( 102, file=file_name )
    do i_p0 = 1, n_p0
        write( 102, '(f15.12,x,f15.12)' ), px(i_p0), pz(i_p0)
    end do
    close( 102 )

    if( rank .eq. 0 ) then
        call pulse_plot
    end if
    
    ! estimate number of trajectories
    allocate( ts_guess( n_p0, LMS_MAX_COUNT ) )
    allocate( n_ts( n_p0 ) )
    do i_p0 = 1, n_p0
        call set_p0( px(i_p0), pz(i_p0) )

! here is only for the first saddle point
        ! sort according to time
        call local_minima_x( LMS_RE_LOWER, LMS_RE_UPPER, LMS_IM_LOWER, LMS_IM_UPPER, SPE, &
              n_ts( i_p0 ), ts_guess( i_p0, : ) )
!!$        ! sort according to importance
!!$        call local_minima_y( LMS_RE_LOWER, LMS_RE_UPPER, LMS_IM_LOWER, LMS_IM_UPPER, SPE, &
!!$              n_ts( i_p0 ), ts_guess( i_p0, : ) )
    end do

    call write_time( 101 )
    write( 101, '(i,10x,a)' ), n_p0, '# n_p0'
    write( 101, '(i,10x,a)' ), p0_start, '# p0_start_index'
    write( 101, '(i,10x,a)' ), p0_start + n_p0 - 1, '# p0_end_index'
    write( 101, '(i,10x,a)' ), sum(n_ts), '# n_traj'
    write( 101, '(a)'), ''
    do i_p0 = p0_start, p0_start+n_p0-1
        index = i_p0 - p0_start + 1
        write( 101, '(i,2x,i)' ), i_p0, n_ts(index)
    end do

    ! start calculation for each trajectory and record of information
    write( file_name, '(a,i3,a)' ), 'dat/traj_', rank, '.dat'
    open( 103, file=file_name )
    do i_p0 = 1, n_p0
        do i_ts = 1, n_ts( i_p0 )
        !do i_ts = 2, 3
            call propagate_with_single_p0( px(i_p0), pz(i_p0), &
                 ts_guess(i_p0, i_ts), ts, Mp, x0, z0, px_inf, pz_inf, &
                 L, n_step, err_spe, ierr, tag )
            write(103, '(2(f0.7,1x), 2(f0.3,1x), 2(f0.2,1x), 2(f0.4,1x), 3(es10.3,1x), i0,1x,i0)' ), &
                 px(i_p0), pz(i_p0), ts, x0, z0, px_inf, pz_inf, L, Mp, n_step, ierr
        end do
    end do
    close( 103 )

    ! finalize
    deallocate( px )
    deallocate( pz )
    deallocate( ts_guess )
    deallocate( n_ts )

    call write_time( 101 )
    close( 101 )
    return
end subroutine slave_thread

integer function amount_of_job( n_data, n_node, rank )
    implicit none
    integer, intent(in):: n_data, n_node, rank
    integer:: n_job, n_res
    n_job = n_data / n_node
    n_res = mod( n_data, n_node )
    if( n_res .eq. 0 ) then
        amount_of_job = n_job
    else
        if( rank < n_res ) then
            amount_of_job = n_job + 1
        else
            amount_of_job = n_job
        end if
    end if
end function amount_of_job

integer function start_of_job( n_data, n_node, rank ) result(job_start)
    implicit none
    integer, intent(in):: n_data, n_node, rank
    integer:: i_rank
    integer, external:: amount_of_job
    job_start = 0
    do i_rank = 0, rank-1
        job_start = job_start + amount_of_job( n_data, n_node, i_rank )
    end do
end function start_of_job

subroutine write_time( os )
    implicit none
    integer, intent(in):: os
    character(12):: real_clock(3);
    integer:: date_time(8);

    call date_and_time( real_clock(1), real_clock(2), &
          real_clock(3), date_time );
    write( os, '(i4,a,i2,a,i2,4x,i2,a,i2,a,i2)'), &
          date_time(1), '/', date_time(2), '/', date_time(3), &
          date_time(5), ':', date_time(6), ':', date_time(7);
    return
end subroutine write_time
