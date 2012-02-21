program spec
    implicit none
    integer, parameter:: fig_nx = 400
    integer, parameter:: fig_nz = 1600
    integer, parameter:: nx = 400
    integer, parameter:: nz = 1600
    integer, parameter:: n_type = 5
    integer, parameter:: n_rank = 501
    double precision, parameter:: grid_lower_x = 0d0
    double precision, parameter:: grid_upper_x = 1d0
    double precision, parameter:: grid_lower_z = -2d0
    double precision, parameter:: grid_upper_z = 2d0
    double precision, parameter:: d_px = ( grid_upper_x - grid_lower_x ) / nx
    double precision, parameter:: d_pz = ( grid_upper_z - grid_lower_z ) / nz
    integer:: n_traj( n_rank ), time_cost( n_rank )
    integer:: traj_count( nx, nz, n_type )
    integer:: err_count( 4 )
    double precision:: cls_Mp( nx, nz, n_type )
    double complex:: qtm_Mp( nx, nz, n_type )
    integer:: i_rank, i_px, i_pz, i_type, ierr
    character(len=128):: file_name 
    double complex:: grid_Mp( nx, nz )
    character(*), parameter:: dir_dat = "../../dat/"
    
    traj_count = 0
    err_count = 0
    cls_Mp = 0d0
    qtm_Mp = (0d0,0d0)

    do i_rank = 0, n_rank-1
        
        write(*,'(i4,a,i4,a)'), i_rank, ' of ', n_rank-1, ' is being proc ...'
        call read_info( i_rank, n_traj(i_rank), time_cost(i_rank), ierr )
        if( ierr .eq. 1 ) then
            write(*,*), "info error, skip read_traj: ", i_rank
            cycle
        end if
        write(*,'(a,i)'), 'n_traj = ', n_traj(i_rank)
        call read_traj( nx, nz, grid_lower_x, d_px, grid_lower_z, d_pz, &
              n_type, i_rank, n_traj(i_rank), traj_count, err_count, &
              qtm_Mp, cls_Mp )
        write(*,'(a)'), repeat('-', 79)

    end do
    
    ! output parallel related information
    call write_rank_info( n_rank, n_traj, time_cost )

    ! output statistical information
    write(*,'(a,i)'), '#valid = ', sum( traj_count(:,:,:) )
    write(*,'(a,i)'), '#T1 = ', sum( traj_count(:,:,1) )
    write(*,'(a,i)'), '#T2 = ', sum( traj_count(:,:,2) )
    write(*,'(a,i)'), '#T3 = ', sum( traj_count(:,:,3) )
    write(*,'(a,i)'), '#T4 = ', sum( traj_count(:,:,4) )
    write(*,'(a,i)'), '#T5 = ', sum( traj_count(:,:,5) )

    write(*,'(a,i)'), '#ts_err',  err_count(1)
    write(*,'(a,i)'), '#rk4_err', err_count(2)
    write(*,'(a,i)'), '#bound',   err_count(3)
    write(*,'(a,i)'), '#caught',  err_count(4)

    
    ! the quantum spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( qtm_Mp(i_px,i_pz,1:4) );
    file_name = dir_dat // 'spec_qtm_all.dat';
    call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
          grid_Mp, fig_nx, fig_nz, file_name )

    ! the classical spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( cls_Mp(i_px,i_pz,1:4) );
    file_name = dir_dat // 'spec_cls_all.dat';
    call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
          grid_Mp, fig_nx, fig_nz, file_name )

    ! the number distribution including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( traj_count(i_px,i_pz,:) )
    file_name = dir_dat // 'spec_number.dat'
    call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
          grid_Mp, fig_nx, fig_nz, file_name )

    ! the quantum spectra for each type of traj
    do i_type = 1, n_type
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = qtm_Mp(i_px,i_pz,i_type);
        write(file_name, '(a,i2,a)'), dir_dat // 'spec_qtm_', i_type, '.dat';
        call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
              grid_Mp, fig_nx, fig_nz, file_name )
    end do

    ! the classical spectra for each type of traj
    do i_type = 1, n_type
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = cls_Mp(i_px,i_pz,i_type);
        write(file_name, '(a,i2,a)'), dir_dat // 'spec_cls_', i_type, '.dat';
        call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
              grid_Mp, fig_nx, fig_nz, file_name )
    end do

    ! the number distribution for each type of traj
    do i_type = 1, n_type
        forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = traj_count(i_px,i_pz,i_type)
        write(file_name, '(a,i2,a)'), dir_dat // 'number_', i_type, '.dat'
        call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
              grid_Mp, fig_nx, fig_nz, file_name )
    end do

    ! the quantum spectra for T1 + T2
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( qtm_Mp(i_px,i_pz,1:2) );
    file_name = dir_dat // 'spec_qtm_1+2.dat';
    call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
          grid_Mp, fig_nx, fig_nz, file_name )

    ! the quantum spectra for T1 + T2 + T3
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( qtm_Mp(i_px,i_pz,1:3) );
    file_name = dir_dat // 'spec_qtm_1+2+3.dat'
    call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
          grid_Mp, fig_nx, fig_nz, file_name )


    write(*,*), 'Hasta la vista'

end program spec

subroutine write_rank_info( n_rank, n_traj, time_cost )
    implicit none
    integer, parameter:: fid_rank_info = 104
    integer, intent(in):: n_rank
    integer, intent(in):: n_traj(n_rank), time_cost(n_rank)
    integer:: i_rank

    open( fid_rank_info, file='../../dat/rank.info')
    do i_rank = 0, n_rank-1
        write( fid_rank_info, '(i4,x,i10,x,i)' ), i_rank, n_traj(i_rank), time_cost(i_rank)
    end do
    close( fid_rank_info )
    return
end subroutine write_rank_info
