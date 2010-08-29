subroutine manage_and_send_data_dbl( i_arry )

    use mod_mmff_core_para, only: n_resc, i_rank, queu_info, arry_info
    use mod_log, only: write_log;

    implicit none;
    integer, intent(in):: i_arry
    integer:: i_queu, i_proc, tag
    double precision, pointer:: send_data(:)
    character(len=200):: text

    interface

        subroutine create_mpi_data_dbl( n_proc, i_proc, n_assign, mpi_data)
            integer, intent(in):: n_proc, i_proc, n_assign(n_proc);
            double precision, pointer:: mpi_data(:)
        end subroutine create_mpi_data_dbl

        subroutine transfer_mpi_data_dbl( n_proc, i_proc, n_assign, data, mpi_data, flag )
            integer, intent(in):: n_proc, i_proc, flag
            integer, intent(in):: n_assign(n_proc);
            double precision, pointer:: data(:);
            double precision, pointer:: mpi_data(:)
        end subroutine transfer_mpi_data_dbl

        subroutine delete_mpi_data_dbl( mpi_data ) 
            double precision, pointer:: mpi_data(:);
        end subroutine delete_mpi_data_dbl

    end interface



    i_queu = arry_info(i_arry) % queu;

    if ( i_rank == 0 ) then
        
        call write_log( repeat('-', 56) );
        write( text, '(2(a,i6), a)' ), &
              'sending queue', i_queu, '    array', i_arry, &
              '     (' // trim(arry_info(i_arry)%name) // ')'; 
        call write_log( text );

        ! for a given data array, loop processes in use
        do i_proc = 1, n_resc
            
            if( queu_info(i_queu) % a_task(i_proc) == 0 ) cycle;
            
            call generate_tag( arry_info(i_arry) % id, i_proc, tag );

            call create_mpi_data_dbl( n_resc, i_proc, queu_info(i_queu) % a_task, send_data );
            call transfer_mpi_data_dbl( n_resc, i_proc, queu_info(i_queu) % a_task, &
                  arry_info(i_arry) % arry_dbl % ptr, send_data, 1 );
            call send_data_mpi_dbl( i_proc, &
                  queu_info(i_queu) % a_task(i_proc), &
                  send_data, &
                  tag, 1 );
            call delete_mpi_data_dbl( send_data );

        end do

    else
        ! for slave, there is no loop; the data array 
        ! can be directly send/recv without splitting

        call generate_tag( arry_info(i_arry) % id, i_rank, tag );

        call send_data_mpi_dbl( 0, &
              queu_info(i_queu) % n_data, &
              arry_info(i_arry) % arry_dbl % ptr, &
              tag, 0 );

    end if

    return;

end subroutine manage_and_send_data_dbl





subroutine manage_and_recv_data_dbl( i_arry )

    use mod_mmff_core_para, only: n_resc, i_rank, queu_info, arry_info
    use mod_log, only: write_log;


    implicit none;
    integer, intent(in):: i_arry
    integer:: i_proc, i_queu, tag
    double precision, pointer:: recv_data(:);
    character(len=200):: text

    interface

        subroutine create_mpi_data_dbl( n_proc, i_proc, n_assign, mpi_data)
            integer, intent(in):: n_proc, i_proc, n_assign(n_proc);
            double precision, pointer:: mpi_data(:)
        end subroutine create_mpi_data_dbl

        subroutine transfer_mpi_data_dbl( n_proc, i_proc, n_assign, data, mpi_data, flag )
            integer, intent(in):: n_proc, i_proc, flag
            integer, intent(in):: n_assign(n_proc);
            double precision, pointer:: data(:);
            double precision, pointer:: mpi_data(:)
        end subroutine transfer_mpi_data_dbl

        subroutine delete_mpi_data_dbl( mpi_data ) 
            double precision, pointer:: mpi_data(:);
        end subroutine delete_mpi_data_dbl

    end interface


    i_queu = arry_info(i_arry) % queu;

    if( i_rank == 0 ) then

        call write_log( repeat('-', 56) );
        write( text, '(2(a,i6), a)' ), &
              'recving queue', i_queu, '    array', i_arry, &
              '     (' // trim(arry_info(i_arry)%name) // ')';  
        call write_log( text );

        ! for a given data array, loop processes in use
        do i_proc = 1, n_resc

            if( queu_info(i_queu) % a_task(i_proc) == 0 ) cycle;

            call generate_tag( arry_info(i_arry) % id, i_proc, tag );

            call create_mpi_data_dbl( n_resc, i_proc, queu_info(i_queu) % a_task, recv_data );
            call recv_data_mpi_dbl( i_proc, &
                  queu_info(i_queu) % a_task(i_proc), &
                  recv_data, &
                  tag, 1 );
            call transfer_mpi_data_dbl( n_resc, i_proc, queu_info(i_queu) % a_task, &
                  arry_info(i_arry) % arry_dbl % ptr, recv_data, -1 );

            call delete_mpi_data_dbl( recv_data );

        end do


    else
        call generate_tag( arry_info(i_arry) % id, i_rank, tag );
        call recv_data_mpi_dbl( 0, &
              queu_info(i_queu) % n_data, &
              arry_info(i_arry) % arry_dbl % ptr, &
              tag, 0 );

    end if

end subroutine manage_and_recv_data_dbl
