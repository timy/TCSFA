subroutine read_info( rank, n_traj, time_cost, ierr )
    implicit none
    integer, intent(in):: rank
    character(len=64):: file_name
    character(len=256):: string
    integer:: n_traj
    integer:: n_p0, p0_start, p0_end
    integer:: time_0(6)
    integer:: time(6)
    integer:: time_diff(6), time_cost, ierr, ierr_time, i
    write( file_name, '(a,i3,a)' ), '../../../dat/info_', rank, '.dat'
    !write( *, '(a)' ), 'reading info from ' // trim( file_name )
    open( 101, file=file_name )
    call read_time( 101, time_0, ierr_time )

    ierr = 0
    read( 101, '(a)', end=201 ), string ! read time 2
    read( 101, '(i,10x,a)', end=201 ), n_p0, string
    read( 101, '(i,10x,a)', end=201 ), p0_start, string
    read( 101, '(i,10x,a)', end=201 ), p0_end, string
    read( 101, '(i,10x,a)', end=201 ), n_traj, string
    read( 101, '(a)', end=201 ), string ! blank line
    do i = 1, n_p0
        read( 101, '(a)', end=201 ), string ! #traj for each p
    end do
    call read_time( 101, time, ierr_time )
    close( 101 )

    !print*, time_0
    !print*, time
    if( ierr_time .eq. 0 ) then
        call time_difference( time_0, time, time_diff )
        time_cost = 12*time_diff(1) + time_diff(2)
        time_cost = 30*time_cost + time_diff(3)
        time_cost = 24*time_cost + time_diff(4)
        time_cost = 60*time_cost + time_diff(5)
        time_cost = 60*time_cost + time_diff(6)
    else
        time_cost = -1
    end if

    return

201 write(*, *), 'the info_file not yet finished ', rank
    ierr = 1
    close( 101 )

    return


end subroutine read_info



subroutine read_time( os, time, ierr_time )
    implicit none
    integer, intent(in):: os
    integer, intent(out):: time(6)
    character(len=1):: s(4)
    integer:: ierr_time
    ierr_time = 0
    read( os, '(i4,a,i2,a,i2,4x,i2,a,i2,a,i2)', end=110 ), &
          time(1), s(1), time(2), s(2), time(3), time(4), &
          s(3), time(5), s(4), time(6)
    return
110 ierr_time = 1 ! if there is no info about the time here
    return
end subroutine read_time




subroutine time_difference( time_start, time_end, time_diff )
    
    implicit none;
    integer, intent(in):: time_start(6), time_end(6)
    integer, intent(out):: time_diff(6)
    integer:: i

    forall(i=1:6) time_diff(i) = time_end(i) - time_start(i);

    if ( time_diff(6) < 0 ) then           ! second
        time_diff(6) = 60 + time_diff(6);
        time_diff(5) = time_diff(5) - 1;
    end if
    if ( time_diff(5) < 0 ) then           ! minute
        time_diff(5) = 60 + time_diff(5);
        time_diff(4) = time_diff(4) - 1;
    end if
    if ( time_diff(4) < 0 ) then           ! hour
        time_diff(4) = 24 + time_diff(4);
        time_diff(3) = time_diff(3) - 1;
    end if
    if ( time_diff(3) < 0 ) then           ! day
        select case ( time_start(2) )
            case(1, 3, 5, 7, 8, 10, 12)
                time_diff(3) = 30 + time_diff(3);
            case(4, 6, 9, 11)
                time_diff(3) = 29 + time_diff(3);
            case(2)
                time_diff(3) = 28 + time_diff(3);
        end select
        time_diff(2) = time_diff(2) - 1;
    end if
    if ( time_diff(2) < 0 ) then           ! month
        time_diff(2) = 12 + time_diff(2);
        time_diff(1) = time_diff(1) - 1;
    end if
    if ( time_diff(1) < 0 ) then
        stop "What the hell happened to the universe?";
    end if

    return;
    
end subroutine time_difference
