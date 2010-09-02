program main

    implicit none;

    integer, parameter:: n_px = 300
    integer, parameter:: n_pz = 600
    integer, parameter:: n_data = n_px*n_pz
    double precision:: grid_px(n_data), grid_pz(n_data), grid_w(n_data)
    integer:: i_px, i_pz
    character(len=200):: fmt


    open(101, file='../data/spec_qtm_all.dat');
    do i_px = 1, n_px
        do i_pz = 1, n_pz
            read( 101, '(3(e15.8,2x))' ), grid_pz(i_pz), grid_px(i_px), &
                  grid_w((i_px-1)*n_pz+i_pz);            
        end do
    end do
    close(101)


    open(101, file='../data/spec_idl.dat');
    write(fmt,'(a,i4,a)'), '(', n_pz, '(e15.8,1x))';
    do i_px = n_px, 1, -1
!        write(101, trim(fmt)), (log10(grid_w((i_px-1)*n_pz+i_pz)), i_pz=1,n_pz);
        write(101, trim(fmt)), ((grid_w((i_px-1)*n_pz+i_pz)), i_pz=1,n_pz);
    end do
    close(101);
    
end program main

