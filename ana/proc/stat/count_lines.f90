subroutine count_lines( fid, n_line )
    implicit none
    integer:: fid, n_line

    n_line = 0
    do
        read( fid, '(a)', end=101 )
        n_line = n_line + 1
    end do

101 rewind( fid )
    return
end subroutine count_lines
