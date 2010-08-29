! MMFF: MPI Framework by Fortran


! global initialization
subroutine mmff_init( mpi_info )
    use mod_mmff_mpi_info, only: T_mmff_mpi_info;
    use mod_mmff_core_para, only:  arry_count, queu_count, n_resc, i_rank

    implicit none;
    type(T_mmff_mpi_info):: mpi_info;
    integer:: i;

    arry_count = 0;
    queu_count = 0;
    n_resc = mpi_info % n_proc - 1;
    i_rank = mpi_info % i_rank;

    return;
end subroutine mmff_init



subroutine mmff_final( )
    use mod_mmff_core_para, only: queu_count, queu_info;

    implicit none;
    integer:: i;

    do i = 1, queu_count
        deallocate( queu_info(i) % a_task );
    end do

    return;
end subroutine mmff_final


! initialize a queue (several arrays with the same size compose a queue)
subroutine mmff_reset( n_data, b_verbose, h_queu )
    use mod_mmff_core_para, only: queu_count, arry_count, i_rank, n_resc, queu_info;

    implicit none;
    integer:: n_data, h_queu
    integer, intent(in):: b_verbose
    integer:: i;

    queu_count = queu_count + 1;

    allocate( queu_info(queu_count) % a_task( n_resc ) );

    if ( b_verbose == 1 ) then
        call assign_workload_with_log( n_data, n_resc, &
              queu_info(queu_count) % a_task );
    else
        call assign_workload( n_data, n_resc, &
              queu_info(queu_count) % a_task );
    end if

    queu_info(queu_count) % id = queu_count;
    queu_info(queu_count) % i_arry = arry_count + 1;
    queu_info(queu_count) % b_arry = 0;
    if( i_rank == 0 ) then
        queu_info(queu_count) % n_data = n_data;
    else
        queu_info(queu_count) % n_data = &
              queu_info(queu_count) % a_task(i_rank);
    end if
    
    n_data = queu_info(queu_count) % n_data;
    h_queu = queu_info(queu_count) % id;

    return;
end subroutine mmff_reset






subroutine mmff_send_data( h_queu )
    
    use mod_mmff_core_para, only: arry_info, queu_info, queu_count, arry_count, &
          MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_INTEGER;

    implicit none;
    integer, intent(in):: h_queu
    integer:: i, i_arry_start, i_arry_end

    if( queu_info(h_queu) % i_arry == 0) return; ! no data array in the queue
    
    i_arry_start = queu_info(h_queu) % i_arry;
    if( h_queu == queu_count ) then
        i_arry_end = arry_count;
    else
        i_arry_end = queu_info(h_queu+1) % i_arry - 1;
    end if

    ! loop all the data arrays in a specified queue
    do i = i_arry_start, i_arry_end

        select case( arry_info(i) % type )
        case( MPI_DOUBLE_PRECISION )
            call manage_and_send_data_dbl( i );
        case( MPI_DOUBLE_COMPLEX )
            call manage_and_send_data_dcp( i );
        case( MPI_INTEGER )
            call manage_and_send_data_int( i );
        case default
        end select

    end do

    return;

end subroutine mmff_send_data





subroutine mmff_recv_data( h_queu )

    use mod_mmff_core_para, only: arry_info, queu_info, queu_count, arry_count, &
          MPI_DOUBLE_PRECISION, MPI_DOUBLE_COMPLEX, MPI_INTEGER;

    implicit none;
    integer, intent(in):: h_queu
    integer:: i, i_arry_start, i_arry_end

    if( queu_info(h_queu) % i_arry == 0) return; ! no data array in the queue
    
    i_arry_start = queu_info(h_queu) % i_arry;
    if( h_queu == queu_count ) then
        i_arry_end = arry_count;
    else
        i_arry_end = queu_info(h_queu+1) % i_arry - 1;
    end if


    do i = i_arry_start, i_arry_end

        select case( arry_info(i) % type )
        case( MPI_DOUBLE_PRECISION )
            call manage_and_recv_data_dbl( i );
        case( MPI_DOUBLE_COMPLEX )
            call manage_and_recv_data_dcp( i );
        case( MPI_INTEGER )
            call manage_and_recv_data_int( i );
        case default
        end select

    end do

    return;

end subroutine mmff_recv_data







subroutine mmff_create_data_dbl( a, f, str_a, h_arry )
    use mod_mmff_core_para, only: queu_count, arry_count, arry_info, queu_info, MPI_DOUBLE_PRECISION

    implicit none;
    double precision, pointer:: a(:);
    external:: f;
    integer:: i, h_arry;
    character(len=*):: str_a;


    arry_count = arry_count + 1;
    arry_info(arry_count) % id = arry_count;
    arry_info(arry_count) % type = MPI_DOUBLE_PRECISION;
    arry_info(arry_count) % queu = queu_count;
    arry_info(arry_count) % name = str_a;

    allocate( a( queu_info(queu_count) % n_data ) );

    do i = 1, queu_info(queu_count) % n_data
        call f( i, a(i) );
    end do

    if( associated(arry_info(arry_count) % arry_dbl % ptr) ) then
        nullify( arry_info(arry_count) % arry_dbl % ptr ); 
    end if
    arry_info(arry_count) % arry_dbl % ptr => a;

    arry_info(arry_count) % b_stat = 1;
    queu_info(queu_count) % b_arry = 1;

    h_arry = arry_info(arry_count) % id;

    return;
end subroutine mmff_create_data_dbl




subroutine mmff_delete_data_dbl( a, h_arry )
    use mod_mmff_core_para, only: arry_info

    implicit none;
    integer, intent(in):: h_arry;
    double precision, pointer:: a(:);

    ! deattach the pointer
    if( associated( arry_info(h_arry) % arry_dbl % ptr ) ) then
        nullify( arry_info(h_arry) % arry_dbl % ptr ); 
    end if
    ! free the heap
    deallocate( a );
    ! toggle the status
    arry_info(h_arry) % b_stat = 0;

    return;

end subroutine mmff_delete_data_dbl








subroutine mmff_create_data_dcp( a, f, str_a, h_arry )
    use mod_mmff_core_para, only: queu_count, arry_count, arry_info, queu_info, MPI_DOUBLE_COMPLEX

    implicit none;
    double complex, pointer:: a(:);
    external:: f;
    integer:: i, h_arry;
    character(len=*):: str_a;


    arry_count = arry_count + 1;
    arry_info(arry_count) % id = arry_count;
    arry_info(arry_count) % type = MPI_DOUBLE_COMPLEX;
    arry_info(arry_count) % queu = queu_count;
    arry_info(arry_count) % name = str_a;

    allocate( a( queu_info(queu_count) % n_data ) );

    do i = 1, queu_info(queu_count) % n_data
        call f( i, a(i) );
    end do

    if( associated(arry_info(arry_count) % arry_dcp % ptr) ) then
        nullify( arry_info(arry_count) % arry_dcp % ptr ); 
    end if
    arry_info(arry_count) % arry_dcp % ptr => a;

    arry_info(arry_count) % b_stat = 1;
    queu_info(queu_count) % b_arry = 1;

    h_arry = arry_info(arry_count) % id;

    return;
end subroutine mmff_create_data_dcp




subroutine mmff_delete_data_dcp( a, h_arry )
    use mod_mmff_core_para, only: arry_info

    implicit none;
    integer, intent(in):: h_arry;
    double complex, pointer:: a(:);

    ! deattach the pointer
    if( associated( arry_info(h_arry) % arry_dcp % ptr ) ) then
        nullify( arry_info(h_arry) % arry_dcp % ptr ); 
    end if
    ! free the heap
    deallocate( a );
    ! toggle the status
    arry_info(h_arry) % b_stat = 0;

    return;

end subroutine mmff_delete_data_dcp









subroutine mmff_create_data_int( a, f, str_a, h_arry )
    use mod_mmff_core_para, only: queu_count, arry_count, arry_info, queu_info, MPI_INTEGER

    implicit none;
    integer, pointer:: a(:);
    external:: f;
    integer:: i, h_arry;
    character(len=*):: str_a;


    arry_count = arry_count + 1;
    arry_info(arry_count) % id = arry_count;
    arry_info(arry_count) % type = MPI_INTEGER;
    arry_info(arry_count) % queu = queu_count;
    arry_info(arry_count) % name = str_a;

    allocate( a( queu_info(queu_count) % n_data ) );

    do i = 1, queu_info(queu_count) % n_data
        call f( i, a(i) );
    end do

    if( associated(arry_info(arry_count) % arry_int % ptr) ) then
        nullify( arry_info(arry_count) % arry_int % ptr ); 
    end if
    arry_info(arry_count) % arry_int % ptr => a;

    arry_info(arry_count) % b_stat = 1;
    queu_info(queu_count) % b_arry = 1;

    h_arry = arry_info(arry_count) % id;

    return;
end subroutine mmff_create_data_int




subroutine mmff_delete_data_int( a, h_arry )
    use mod_mmff_core_para, only: arry_info

    implicit none;
    integer, intent(in):: h_arry;
    integer, pointer:: a(:);

    ! deattach the pointer
    if( associated( arry_info(h_arry) % arry_int % ptr ) ) then
        nullify( arry_info(h_arry) % arry_int % ptr ); 
    end if
    ! free the heap
    deallocate( a );
    ! toggle the status
    arry_info(h_arry) % b_stat = 0;

    return;

end subroutine mmff_delete_data_int




subroutine mmff_get_task( h_queu, a_task )
    use mod_mmff_core_para, only: n_resc, queu_info

    implicit none;
    integer, intent(in):: h_queu
    integer, pointer:: a_task(:)
    integer:: i

    allocate( a_task( n_resc ) );
    forall(i=1:n_resc) a_task(i) = queu_info(h_queu) % a_task(i);

    return;
end subroutine mmff_get_task



subroutine mmff_set_task( h_queu, a_task )
    use mod_mmff_core_para, only: queu_count, arry_count, i_rank, n_resc, queu_info;

    implicit none;
    integer, intent(in):: h_queu
    integer, pointer:: a_task(:)
    integer:: i;

    if( sum( a_task ) == sum( queu_info(h_queu) % a_task ) ) then
        forall(i=1:n_resc) queu_info(h_queu) % a_task(i) = a_task(i);
    else
        stop 'invalid assign in mmff_set_task'
    end if

    if( i_rank > 0 ) then
        queu_info(h_queu) % n_data = queu_info(h_queu) % a_task(i_rank);
    end if

    deallocate( a_task );

    return;
end subroutine mmff_set_task




subroutine mmff_broadcast_master_int( na, a )
    use mod_mmff_core_para, only: n_resc

    implicit none
    integer, intent(in):: na, a(na)
    integer:: first, last, n_data, i_proc
    integer:: h_q1, h_q2, h_a1, h_a2
    integer, external:: null_int
    integer, pointer:: data_size(:), data(:)
    
    interface
        subroutine mmff_create_data_int( a, f, str, h_arry )
            integer, pointer:: a(:)
            external:: f;
            character(len=*):: str;
            integer, intent(in):: h_arry;
        end subroutine mmff_create_data_int

        subroutine mmff_delete_data_int( a, h_arry )
            integer, pointer:: a(:)
            integer, intent(in):: h_arry
        end subroutine mmff_delete_data_int
    end interface


    ! first, send size of array to each proc
    n_data = n_resc;
    call mmff_reset( n_data, 0, h_q1 );
    call mmff_create_data_int( data_size, null_int, 'broadcast_size', h_a1 );
    forall(i_proc=1:n_resc) data_size(i_proc) = na;
    call mmff_send_data( h_q1 );
    call mmff_delete_data_int( data_size, h_a1 );

    !second, send the actual data array to each proc
    n_data = n_resc * na;
    call mmff_reset( n_data, 0, h_q2 );
    call mmff_create_data_int( data, null_int, 'broadcast_data', h_a2 );
    do i_proc = 1, n_resc
        first = (i_proc - 1) * na + 1;
        last = first + na - 1;
        data(first:last) = a(1:na);
    end do
    call mmff_send_data( h_q2 );
    call mmff_delete_data_int( data, h_a2 );

    return;
end subroutine mmff_broadcast_master_int





subroutine mmff_broadcast_slave_int( a )
    use mod_mmff_core_para, only: n_resc

    implicit none
    integer:: n_data, na, i
    integer:: h_q1, h_q2, h_a1, h_a2
    integer, external:: null_int
    integer, pointer:: data_size(:), data(:), a(:)
    
    interface
        subroutine mmff_create_data_int( a, f, str, h_arry )
            integer, pointer:: a(:)
            external:: f;
            character(len=*):: str;
            integer, intent(in):: h_arry;
        end subroutine mmff_create_data_int

        subroutine mmff_delete_data_int( a, h_arry )
            integer, pointer:: a(:)
            integer, intent(in):: h_arry
        end subroutine mmff_delete_data_int
    end interface

    ! first, recv size of array in each proc
    n_data = n_resc;
    call mmff_reset( n_data, 0, h_q1 );
    call mmff_create_data_int( data_size, null_int, 'braodcast_size', h_a1 );
    call mmff_recv_data( h_q1 );
    na = data_size(1);
    allocate( a(na) );
    call mmff_delete_data_int( data_size, h_a1 );

    ! second, recv the actual data array to each proc
    n_data = na * n_resc;
    call mmff_reset( n_data, 0, h_q2 );
    call mmff_create_data_int( data, null_int, 'broadcast_data', h_a2 );
    call mmff_recv_data( h_q2 );
    forall(i=1:n_data) a(i) = data(i);
    call mmff_delete_data_int( data, h_a2 );

    return;
end subroutine mmff_broadcast_slave_int




subroutine null_dbl( x, y )
    implicit none;
    integer, intent(in):: x
    double precision:: y

    y = 0d0;

    return;
end subroutine null_dbl

subroutine null_dcp( x, y )
    implicit none;
    integer, intent(in):: x
    double complex:: y

    y = (0d0, 0d0);

    return;
end subroutine null_dcp

subroutine null_int( x, y )
    implicit none;
    integer, intent(in):: x
    integer:: y

    y = 0;

    return;
end subroutine null_int
