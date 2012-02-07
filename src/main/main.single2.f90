#include '../include/inc_field.h'
#include '../include/inc_ts_guess.h'
#include '../include/inc_grid_boundary.h'
#include '../include/inc_grid_size.h'
#include '../include/inc_rk4.h'

program test
    use mod_pulse
    implicit none

    character(len=64):: file_name

    integer, external:: amount_of_job, start_of_job
    integer:: n_p0, p0_start
    integer:: ix, iz, i_p0, i_ts, index
    double precision, allocatable:: px(:), pz(:)
    integer, allocatable:: n_ts(:)
    double complex, allocatable:: ts_guess(:,:)
    double precision:: rand_px_0, rand_pz_0

    double complex:: ts, Mp
    double precision:: x0, z0, px_inf, pz_inf, L
    integer:: n_pass_x, n_pass_z, ierr, n_near_core
    integer:: seed(1)
    integer, parameter:: tag = 0

    double complex, external:: SPE

    ! organize my tasks
    n_p0 = N_P0
    p0_start = 1

    ! init the momentum according to my rank
    allocate( px(n_p0) )
    allocate( pz(n_p0) )
    seed(1) = 517
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


!        px(index) = 0.375108648503d0    
!        pz(index) = 0.115438116627d0

    end do

    ! output px, pz
    write( file_name, '(a)' ), 'dat/p.dat'
    open( 102, file=file_name )
    do i_p0 = 1, n_p0
        write( 102, '(f15.12,x,f15.12)' ), px(i_p0), pz(i_p0)
    end do
    close( 102 )

    ! set laser pulse
    call set_pulse( E0, OM, NC, XI, PH )
    call set_pulse_t0( (1d0, 0d0) )
    
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
    write( file_name, '(a)' ), 'dat/traj.dat'
    open( 103, file=file_name )
    do i_p0 = 1, n_p0
        do i_ts = 1, n_ts( i_p0 )
        !do i_ts = 2, 3
        
            
            write(*,'(4(e15.7,1x))'), px(i_p0), pz(i_p0), ts 
            call propagate_with_single_p0( px(i_p0), pz(i_p0), &
                  ts_guess(i_p0, i_ts), ts, Mp, x0, z0, px_inf, pz_inf, &
                  L, n_pass_x, n_pass_z, n_near_core, ierr, tag )
            write( 103, '(11(e16.8, 1x), 2(i2,1x))' ), &
                  px(i_p0), pz(i_p0), ts, x0, z0, px_inf, pz_inf, L, &
                  Mp, n_near_core, ierr
            
        end do
    end do
    close( 103 )

    ! finalize
    deallocate( px )
    deallocate( pz )
    deallocate( ts_guess )
    deallocate( n_ts )

end program test

