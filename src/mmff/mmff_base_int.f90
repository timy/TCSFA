subroutine send_data_mpi_int( i_proc, n_elem, send_data, tag, b_verbose )

    use mod_log, only: write_log

    implicit none;
    include 'mpif.h'

    integer, intent(in):: i_proc, n_elem, tag, b_verbose;
    integer:: send_data(n_elem);
    integer:: mpi_ierr;
    character(len=200):: text;

    call mpi_send( send_data, n_elem, MPI_INTEGER, &
          i_proc, tag, MPI_COMM_WORLD, mpi_ierr );

    if( b_verbose ) then

        if ( mpi_ierr == MPI_SUCCESS ) then

            write( text, '(3(a,i6))' ), 'integer data', tag, ' -> node', i_proc, ' : size =', n_elem; 
            call write_log( text );

        else 

            write( text, '(a, i8, a, i8, a, i8)'), &
                  'error from data ', tag, ' in process ', i_proc, &
                  ': ', mpi_ierr;
            call write_log( text );

        end if

    end if

    return;

end subroutine send_data_mpi_int



subroutine recv_data_mpi_int( i_proc, n_elem, recv_data, tag, b_verbose )

    use mod_log, only: write_log

    implicit none;
    include 'mpif.h'

    integer, intent(in):: i_proc, n_elem, tag, b_verbose;
    integer:: recv_data(n_elem);
    integer:: mpi_ierr, status(4);
    character(len=200):: text;

    call mpi_recv( recv_data, n_elem, MPI_INTEGER, i_proc, &
          tag, MPI_COMM_WORLD, status, mpi_ierr );

    if( b_verbose ) then

        if ( mpi_ierr == MPI_SUCCESS ) then

            write( text, '(3(a,i6))' ), 'integer data', tag, ' <- node', i_proc, ' : size =', n_elem; 
            call write_log( text );

        else

            write( text, '(a, i8, a, i8, a, i8)'), &
                  'error from data ', tag, ' in process ', i_proc, &
                  ': ', mpi_ierr;
            call write_log( text );

        end if

    end if
    return;

end subroutine recv_data_mpi_int





! create array 'mpi_data' which is used for sending data to slave process
subroutine create_mpi_data_int( n_proc, i_proc, n_assign, mpi_data )

    implicit none;

    integer, intent(in):: n_proc, i_proc
    integer, intent(in):: n_assign(n_proc)
    integer, pointer:: mpi_data(:)
    integer:: i, start, n_elem

    n_elem = n_assign(i_proc);
    if( associated( mpi_data ) ) deallocate( mpi_data );
    allocate( mpi_data(n_elem) );

    return;

end subroutine create_mpi_data_int



subroutine transfer_mpi_data_int( n_proc, i_proc, n_assign, data, mpi_data, flag )

    implicit none;

    integer, intent(in):: n_proc, i_proc, flag
    integer, intent(in):: n_assign(n_proc)
    integer, pointer:: data(:)
    integer, pointer:: mpi_data(:)
    integer:: i, n_elem, start

    n_elem = n_assign(i_proc);


    if( i_proc == 1 ) then
        start = 0;
    else
        start = sum( n_assign(1:i_proc-1) );
    end if

    if( associated( mpi_data ) ) then

        select case (flag)
        case(1)
            forall(i=1:n_elem) mpi_data(i) = data(start+i);
        case(-1)
            forall(i=1:n_elem) data(start+i) = mpi_data(i);
        case default
        end select

    else
        stop "transfer_mpi_data_int: array not yet allocated!";
    end if

    return;

end subroutine transfer_mpi_data_int






! delete array 'mpi_data'
subroutine delete_mpi_data_int( mpi_data )

    implicit none;

    integer, pointer:: mpi_data(:);

    if( associated( mpi_data ) ) deallocate( mpi_data );

    return;

end subroutine delete_mpi_data_int


