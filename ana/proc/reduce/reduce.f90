program selc
    implicit none
    integer:: selc_count, selc_count_total
    integer, parameter:: n_rank = 701
    integer:: n_traj, time_cost(n_rank), i_rank
    character(*), parameter:: dir_dat = "../../dat/"
    integer, parameter:: fid_output = 131

    selc_count_total = 0
    
    open( fid_output, file='../../dat/reduced_data.dat' )

    do i_rank = 0, n_rank-1

        write(*,'(i4,a,i4,a)'), i_rank, ' of ', n_rank-1, ' is being proc ...'
        call read_info( i_rank, n_traj, time_cost(i_rank) )
        write(*,'(a,i)'), 'n_traj = ', n_traj
        call read_traj( i_rank, n_traj, fid_output, selc_count )
        write(*,'(a,i)'), 'number of selected traj: ', selc_count
        selc_count_total = selc_count_total + selc_count

    end do

    close( fid_output )
    write(*,'(a,i)'), 'reduced number of traj: ', selc_count_total
    write(*,'(a,i)'), 'hasta la vista'

end program selc
