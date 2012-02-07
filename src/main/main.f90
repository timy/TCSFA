program main
    implicit none
    include 'mpif.h'

    integer:: ierr, n_node, rank

    call mpi_init( ierr )
    call mpi_comm_size( MPI_COMM_WORLD, n_node, ierr )
    call mpi_comm_rank( MPI_COMM_WORLD, rank, ierr )
    call slave_thread( n_node, rank )
    call mpi_finalize( ierr )

end program main
