program main
    implicit none
    include 'mpif.h'
    integer, parameter:: fig_nx = 400
    integer, parameter:: fig_nz = 1600
    integer, parameter:: nx = 400
    integer, parameter:: nz = 1600
    integer, parameter:: n_type = 5
    integer, parameter:: n_rank_raw_data = 611
    character(*), parameter:: dir_dat = "../../dat/"
    double precision, parameter:: grid_lower_x = 0d0
    double precision, parameter:: grid_upper_x = 1d0
    double precision, parameter:: grid_lower_z = -2d0
    double precision, parameter:: grid_upper_z = 2d0
    double precision, parameter:: d_px = ( grid_upper_x - grid_lower_x ) / nx
    double precision, parameter:: d_pz = ( grid_upper_z - grid_lower_z ) / nz
    integer:: ierr, n_node, rank, i_rank, index, status(4)
    integer:: n_file, file_start, i_file
    integer, external:: amount_of_job, start_of_job
    integer:: n_traj, time_cost
    integer:: traj_count( nx, nz, n_type ), err_count(4), &
          traj_count_buf( nx, nz, n_type ), err_count_buf(4)
    double complex:: qtm_Mp( nx, nz, n_type ), qtm_Mp_buf( nx, nz, n_type )
    double precision:: cls_Mp( nx, nz, n_type ), cls_Mp_buf( nx, nz, n_type )
    integer, parameter:: fid_info = 101

    ! obtain info n_node and current rank
    call mpi_init( ierr )
    call mpi_comm_size( MPI_COMM_WORLD, n_node, ierr )
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )

    ! organize my tasks
    n_file = amount_of_job( n_rank_raw_data, n_node, rank )
    file_start = start_of_job( n_rank_raw_data, n_node, rank )
    print*, 'rank', rank, 'n_file', n_file, 'start', file_start

    ! initialize the data arrays
    traj_count = 0
    err_count = 0
    cls_Mp = 0d0
    qtm_Mp = dcmplx(0d0,0d0)

    if( rank .eq. 0 ) then
        open( fid_info, file=dir_dat//"spec.info")
    end if

    ! read trajectories from files
    do i_file = 1, n_file
        index = file_start + i_file - 1
        call read_info( index, n_traj, time_cost )
        call read_traj( nx, nz, grid_lower_x, d_px, grid_lower_z, d_pz, &
              n_type, index, n_traj, traj_count, err_count, qtm_Mp, cls_Mp )
        if( rank .eq. 0 ) then
            write( fid_info, '(a,i4,a,i4)'), 'rank 0:', i_file, 'of', n_file 
        end if
    end do

    ! syncronize it 
    call mpi_barrier( MPI_COMM_WORLD )

    ! when all data is ready, send them back to rank 0 and output
    if( rank .eq. 0 ) then
        ! recv data from other ranks, sum them up
        write( fid_info, '(a)' ), 'start to recv data ...'
        do i_rank = 1, n_node

            call mpi_recv( traj_count_buf, nx*nz*n_type, MPI_INTEGER, &
                  i_rank, 1000+i_rank, MPI_COMM_WORLD, status, ierr )
            traj_count = traj_count + traj_count_buf
            call mpi_recv( err_count_buf, 4, MPI_INTEGER, &
                  i_rank, 2000+i_rank, MPI_COMM_WORLD, status, ierr )
            err_count = err_count + err_count_buf
            call mpi_recv( cls_Mp_buf, nx*nz*n_type, MPI_DOUBLE_PRECISION, &
                  i_rank, 3000+i_rank, MPI_COMM_WORLD, status, ierr )
            cls_Mp = cls_Mp + cls_Mp_buf
            call mpi_recv( qtm_Mp_buf, nx*nz*n_type, MPI_DOUBLE_COMPLEX, &
                  i_rank, 4000+i_rank, MPI_COMM_WORLD, status, ierr )
            qtm_Mp = qtm_Mp + qtm_Mp_buf
            
            write( fid_info, '(a,i4,a)'), 'thread ', i_rank, ' done. '

        end do
        ! output results
        call spec_output( nx, nz, n_type, fig_nx, fig_nz, grid_lower_x, grid_lower_z, &
              d_px, d_pz, traj_count, err_count, cls_Mp, qtm_Mp )
        close( fid_info )
    else
        ! send data back to the rank 0
        call mpi_send( traj_count, nx*nz*n_type, MPI_INTEGER, &
              0, 1000+rank, MPI_COMM_WORLD, ierr )
        call mpi_send( err_count, 4, MPI_INTEGER, &
              0, 2000+rank, MPI_COMM_WORLD, ierr )
        call mpi_send( cls_Mp, nx*nz*n_type, MPI_DOUBLE_PRECISION, &
              0, 3000+rank, MPI_COMM_WORLD, ierr )
        call mpi_send( qtm_Mp, nx*nz*n_type, MPI_DOUBLE_COMPLEX, &
              0, 4000+rank, MPI_COMM_WORLD, ierr )
    end if

    call mpi_finalize( ierr )

end program main


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
