program selc
    implicit none
    integer, parameter:: n_sel = 40
    double precision:: z0(n_sel)
    double precision:: x0(n_sel)
    double precision:: r0(n_sel)
    integer:: fid_selc(n_sel), selc_count, selc_count_total
    integer, parameter:: n_rank = 401
    integer:: n_traj(n_rank), time_cost(n_rank), i_rank, i_s
    character(*), parameter:: dir_dat = "../../dat/"
    character(len=64):: filename
    
    call sample_points( z0, x0, r0, n_sel )

    selc_count_total = 0

    do i_s = 1, n_sel
        fid_selc(i_s) = 110 + i_s
        write( filename, '(a,i3,a)' ), dir_dat // 'selc_', i_s, '.dat'
        open( fid_selc(i_s), file=filename )
    end do

    do i_rank = 0, n_rank-1

        write(*,'(i4,a,i4,a)'), i_rank, ' of ', n_rank-1, ' is being proc ...'
        call read_info( i_rank, n_traj(i_rank), time_cost(i_rank) )
        write(*,'(a,i)'), 'n_traj = ', n_traj(i_rank)
        call read_traj( n_sel, x0, z0, r0, i_rank, n_traj(i_rank), fid_selc, selc_count )
        write(*,'(a,i)'), 'number of selected traj: ', selc_count
        selc_count_total = selc_count_total + selc_count

    end do
    
    do i_s = 1, n_sel
        close( fid_selc(i_s) )
    end do

    write(*,'(a,i)'), 'hasta la vista'

end program selc
