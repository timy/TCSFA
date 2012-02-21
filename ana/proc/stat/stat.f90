program stat
    implicit none
    integer, parameter:: i_selc = 6
    integer:: n_traj
    character(*), parameter:: dir_dat = "../../dat/"
    character(len=64):: filename
    integer, parameter::fid_selc = 101

    write( filename, '(a,i3,a)' ), dir_dat // 'selc_', i_selc, '.dat'
    open( fid_selc, file=filename )
    call count_lines( fid_selc, n_traj )
    write(*, '(a,i)'), "n_traj = ", n_traj
    call read_selc( n_traj, fid_selc )
    close( fid_selc )

    write(*, '(a)'), 'done!'

end program stat
