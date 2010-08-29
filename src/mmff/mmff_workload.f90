! ################################################################################
! workload.f90 is used to assign balanced tasks to each process; Two subroutines 
! are included: assign_workload_with_log, assign_workload
! ################################################################################

subroutine assign_workload_with_log( n_data, n_proc, n_assign )

    use mod_log, only: write_log

    implicit none;
    integer, intent(in):: n_data, n_proc;
    integer:: n_assign(n_proc);
    character(len=200):: text;
    integer:: i;

    call assign_workload( n_data, n_proc, n_assign );    

    ! record the information for workload and assignment
    call write_log( repeat( '-', 40 ) );
    write(text, '(a,i9)'),         'amount of all data to go:   ', n_data 
    call write_log(text);
    write(text, '(a,i9)'),         'number of total processes:  ', n_proc;
    call write_log(text);
    write(text, '(a,i9)'),         'number of processes in use:  ', count(n_assign>0);
    call write_log(text);

    do i = 1, n_proc
        write(text, '(a,i5,a,i8)'), 'rank:', i, '->', n_assign(i);
        call write_log(text);
    end do

    return;
end subroutine assign_workload_with_log


subroutine assign_workload( n_data, n_proc, n_assign )

    implicit none;
    integer, intent(in):: n_data, n_proc;
    integer:: n_assign(n_proc)
    integer:: i, nb, nr;

    nb = n_data / n_proc;
    forall(i=1:n_proc) n_assign(i) = nb;
    nr = n_data - n_proc * nb;


    if( nr == n_data ) then
        ! n_data < n_proc
        do i = 1, n_data
            n_assign(i) = 1;
        end do

    elseif( nr == 0) then
        ! n_data == n_proc

    elseif( nr > 0) then
        ! n_data > n_proc && mod(n_data) /= 0
        do i = 1, nr
            n_assign(i) = n_assign(i) + 1;
        end do

    end if

    return;
end subroutine assign_workload

integer function estimate_workload_slave( n_data, n_proc, rank )
    implicit none;
    integer, intent(in):: n_data, n_proc, rank;
    integer:: n_assign(n_proc);
    call assign_workload( n_data, n_proc, n_assign );
    estimate_workload_slave = n_assign(rank);
    return;
end function estimate_workload_slave
