! This module is used to record log, it will automatively add time to the event, 
! and user should provide the information it records.

module mod_log
    implicit none;
    character(*), parameter:: mod_log_filename = 'LOG';
    type, private:: class_log
        character(len=22):: current_time;
        character(len=200):: info
        integer:: opt = 1;
    end type class_log

contains

    subroutine delete_log
        use ifport;
        implicit none;
        logical:: exists;
        integer:: count;

        inquire( file='LOG', exist=exists );
        if ( exists ) count = system('rm LOG');
        return;
    end subroutine delete_log

    subroutine set_log_mode( opt )
        implicit none;
        integer, intent(in):: opt;
        type(class_log):: current_log;

        current_log % opt = opt;
        return;
    end subroutine set_log_mode

    subroutine write_log( info )
        implicit none
        character(*), intent(in):: info;
        type(class_log):: current_log;
        character(12):: real_clock(3);
        integer:: date_time(8);


        open( file=mod_log_filename, unit=99, position='APPEND' );
        if( current_log%opt == 0 ) then
            call date_and_time( real_clock(1), real_clock(2), &
                  real_clock(3), date_time );
            write(current_log%current_time, '(i4,a,i2,a,i2,4x,i2,a,i2,a,i2)'), &
                  date_time(1), '/', date_time(2), '/', date_time(3), &
                  date_time(5), ':', date_time(6), ':', date_time(7);
            write( 99, '(a, 4x, a)'), current_log%current_time, repeat('-', 24);
        end if
        write(current_log%info,'(a)'), info;
        write( 99, '(a)'), current_log%info;
        close( 99 );

        write( *, '(a)'), current_log%info;

        return;
    end subroutine write_log
end module mod_log
