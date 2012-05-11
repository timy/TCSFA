subroutine calc_w( nx, nz, d_px, d_pz, grid_lower_x, grid_lower_z, &
      grid_Mp, fig_nx, fig_nz, file_name )
    implicit none
    integer, parameter:: fid_spec = 103
    integer, intent(in):: nx, nz, fig_nx, fig_nz
    double precision, intent(in):: d_px, d_pz, grid_lower_x, grid_lower_z
    double complex, intent(in):: grid_Mp( nx, nz )
    double precision:: fig_px(fig_nx), fig_pz(fig_nz), fig_w(fig_nx,fig_nz)
    character(len=128):: file_name
    integer:: i_px, i_pz, ni_x, ni_z

    forall(i_px=1:fig_nx) fig_px(i_px) = grid_lower_x + (i_px - 0.5d0) * d_px;
    forall(i_pz=1:fig_nz) fig_pz(i_pz) = grid_lower_z + (i_pz - 0.5d0) * d_pz;
    forall(i_px=1:fig_nx,i_pz=1:fig_nz) fig_w(i_px,i_pz) = 0d0;

    ! # of grids considered for one pixel in the figure
    ni_x = nx / fig_nx; ! nx should be larger than fig_nx
    ni_z = nz / fig_nz; ! nz should be larger than fig_nz
    
    do i_px = 1, fig_nx ! index of x direction in the figure
        do i_pz = 1, fig_nz ! index of z direction in the figure
            
            fig_w(i_px,i_pz) = sum( cdabs( grid_Mp( ((i_px-1)*ni_x+1):(i_px*ni_x), ((i_pz-1)*ni_z+1):(i_pz*ni_z)) )**2 );

        end do
    end do

    ! now write into file
    write(*, '(a, a)'), 'now writing to file: ', file_name;
    open( fid_spec, file=trim(file_name) );
    do i_px = 1, fig_nx
        do i_pz = 1, fig_nz
            
            if ( dabs(fig_w(i_px,i_pz)) < 1d-99 ) fig_w(i_px,i_pz) = 1d-99;
            if ( dabs(fig_w(i_px,i_pz)) > 1d99 ) fig_w(i_px,i_pz) = 1d99;
            write( fid_spec, '(2(f0.4,1x), es10.3)' ), fig_pz(i_pz), fig_px(i_px), fig_w(i_px,i_pz);
                
        end do
    end do
    close( FID_SPEC );

    return;
end subroutine calc_w
