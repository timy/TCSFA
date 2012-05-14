program selc
    implicit none
    include 'sample_number.f'
    double precision:: z0(n_sel)
    double precision:: x0(n_sel)
    double precision:: r0(n_sel)
    integer:: fid_selc(n_sel), selc_count, selc_count_total
    integer, parameter:: n_rank = 601
    integer:: n_traj(n_rank), time_cost(n_rank), i_rank, i_s
    character(*), parameter:: dir_dat = "../../dat/selc/"
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


subroutine sample_points( z0, x0, r0, n )

    implicit none
    integer, intent(in):: n
    double precision:: z0(n), x0(n), r0(n)
    integer:: i

    include 'sample.f'
    write(*, '(a)'), "Trajectories corresponding to the following points will be extracted"
    do i = 0, n
        write(*, '(f15.8, f15.8)'), z0(i), x0(i)
    end do

    ! the radius
    do i = 0, n
        r0(i) = 0.005d0;
    end do

    return
end subroutine sample_points
