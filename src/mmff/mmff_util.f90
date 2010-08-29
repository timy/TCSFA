subroutine generate_tag( arry_id, i_proc, tag )

    implicit none;
    integer, intent(in):: arry_id, i_proc;
    integer:: tag;

    tag = arry_id * 1000 + i_proc

    return;
end subroutine generate_tag



subroutine timer_tic_with_log( time_start )

    use mod_log, only: write_log
    implicit none

    integer:: time_start(6)
    character(len=200):: text;

    call timer_tic( time_start );
    write(text, '(a,6(i4,a))'), 'Starting time is ', &
          time_start(1), 'Y ', time_start(2), 'M ', time_start(3), 'D ', &
          time_start(4), 'h ', time_start(5), 'm ', time_start(6), 's';
    call write_log( text );

    return;

end subroutine timer_tic_with_log




subroutine timer_toc_with_log( time_start )

    use mod_log, only: write_log
    implicit none

    integer, intent(in):: time_start(6)
    integer:: time_diff(6)
    character(len=200):: text;

    call timer_toc( time_start, time_diff );
    write(text, '(a,6(i4,a))'), 'Elapsed time is ', &
          time_diff(1), 'Y ', time_diff(2), 'M ', time_diff(3), 'D ', &
          time_diff(4), 'h ', time_diff(5), 'm ', time_diff(6), 's';
    call write_log( text );

    return;

end subroutine timer_toc_with_log



subroutine timer_tic( time_start )
    
    implicit none;
    integer, intent(out):: time_start(6)

    call get_now( time_start );

    return;

end subroutine timer_tic



subroutine timer_toc( time_start, time_diff )
    
    implicit none;
    integer, intent(in):: time_start(6)
    integer, intent(out):: time_diff(6)
    integer:: time_end(6), i;

    call get_now( time_end );

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
        stop "What the hell happened with the universe?";
    end if

    return;
    
end subroutine timer_toc



subroutine get_now( time_now )

    implicit none;
    integer, intent(out):: time_now(6)
    integer:: month, day, year, now(3)
    integer:: i
    
    call idate( month, day, year );
    call itime( now );

    time_now(1) = year;
    time_now(2) = month;
    time_now(3) = day;
    forall(i=1:3) time_now(i+3) = now(i);

    return;

end subroutine get_now
