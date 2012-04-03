subroutine spec_output( nx, nz, n_type, fig_nx, fig_nz, &
      grid_lower_x, grid_lower_z, d_px, d_pz, traj_count, err_count, cls_Mp, qtm_Mp )
    
    implicit none
    integer, intent(in):: nx, nz, n_type, fig_nx, fig_nz
    double precision, intent(in):: grid_lower_x, grid_lower_z, d_px, d_pz
    integer, intent(in):: traj_count( nx, nz, n_type ), err_count(4)
    double complex, intent(in):: qtm_Mp( nx, nz, n_type )
    double precision, intent(in):: cls_Mp( nx, nz, n_type )
    integer, parameter:: fid_info = 101
    integer, parameter:: fid_spec = 102
    character(*), parameter:: dir_dat = "../../dat/"
    double complex:: grid_Mp( nx, nz )
    integer:: i_px, i_pz, i_type
    character(len=128):: file_name

    open( fid_info, file=dir_dat//'spec.info' )
    ! output statistical information
    write(fid_info,'(a,i)'), '#valid = ', sum( traj_count(:,:,:) )
    write(fid_info,'(a,i)'), '#T1 = ', sum( traj_count(:,:,1) )
    write(fid_info,'(a,i)'), '#T2 = ', sum( traj_count(:,:,2) )
    write(fid_info,'(a,i)'), '#T3 = ', sum( traj_count(:,:,3) )
    write(fid_info,'(a,i)'), '#T4 = ', sum( traj_count(:,:,4) )
    write(fid_info,'(a,i)'), '#T5 = ', sum( traj_count(:,:,5) )
    
    write(fid_info,'(a,i)'), '#ts_err',  err_count(1)
    write(fid_info,'(a,i)'), '#rk4_err', err_count(2)
    write(fid_info,'(a,i)'), '#bound',   err_count(3)
    write(fid_info,'(a,i)'), '#caught',  err_count(4)
    close( fid_info )
        
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
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( traj_count(i_px,i_px,:) )
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
end subroutine spec_output

subroutine spec_map_output( nx, nz, n_type, fig_nx, fig_nz, &
      grid_lower_x, grid_lower_z, d_px, d_pz, err_spe)
    
    implicit none
    integer, intent(in):: nx, nz, n_type, fig_nx, fig_nz
    double precision, intent(in):: grid_lower_x, grid_lower_z, d_px, d_pz
    double precision, intent(in):: err_spe( nx, nz, n_type )
    integer, parameter:: fid_spec = 102
    character(*), parameter:: dir_dat = "../../dat/"
    double complex:: grid_Mp( nx, nz )
    integer:: i_px, i_pz
    character(len=128):: file_name
    
    ! the quantum spectra including all types of traj
    forall(i_px=1:nx,i_pz=1:nz) grid_Mp(i_px,i_pz) = sum( err_spe(i_px,i_pz,1:4) );
    file_name = dir_dat // 'spec_err_spe.dat';
    call calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
          err_spe, fig_nx, fig_nz, file_name )
    
end subroutine spec_map_output
