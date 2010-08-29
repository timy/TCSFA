#define FILE_NAME_SELECT '../dat/select 1.dat'
#define FILE_NAME_FILTER '../dat/filter.dat'
#define FILE_NAME_TRAJ_M '../dat/traj_m.dat'

#define FID_SELECT 121
#define FID_FILTER 122
#define FID_TRAJ_M 123

subroutine read_data_and_plot_traj()

    implicit none
    integer:: io, n_line;
    character(len=200):: line, text
    double precision:: px_0, pz_0, ts_re, ts_im, x_t0, z_t0, &
          px_inf, pz_inf, M_re, M_im
    double precision:: pinf_x, pinf_z, x0, z0
    double complex:: amp_M, ts
    integer:: i_data, ierr, n_pass_x, n_pass_z, i_type

    open ( FID_SELECT, file = FILE_NAME_SELECT, status='OLD' );
    open ( FID_FILTER, file = FILE_NAME_FILTER );
    open ( FID_TRAJ_M, file = FILE_NAME_TRAJ_M );

    n_line = 0;
    do
        read ( FID_SELECT, '(A)', iostat = io) line
        if (io < 0) exit
        if (len_trim(line) == 0) cycle
        n_line = n_line + 1;
    end do
    print *, 'number of trajectories in this bin:', n_line;
    
    rewind( FID_SELECT );
    
    do i_data = 1, n_line

        print*, repeat('-', 40) // 'data', i_data;
        open( 999, file='../dat/temp' );
        write( text, '(a,i3,a)'), '../dat/traj_', i_data, '.dat'
        write( 999, '(a)' ), text;
        close( 999 );

        read(FID_SELECT, '(10(e15.8,1x),4(i2,1x))'),  &
              px_0, pz_0, ts_re, ts_im, &
              x_t0, z_t0, px_inf, pz_inf, &
              M_re, M_im, n_pass_x, n_pass_z, ierr, i_type;

        call propagate_with_single_p0( px_0, pz_0, dcmplx( ts_re, ts_im ), &
              ts, amp_M, x0, z0, pinf_x, pinf_z, n_pass_x, n_pass_z, ierr )

        call set_criteria( i_data, x0, z0, px_0, pz_0, pinf_x, pinf_z, ts, amp_M, ierr, i_type )

!        print*, 'Comparison  - px_inf & pinf_x:', px_inf, pinf_x;
!        print*, 'Comparison  - pz_inf & pinf_z:', pz_inf, pinf_z;

    end do

    
    
    close( FID_FILTER );
    close( FID_SELECT );
    close( FID_TRAJ_M );
    return;
end subroutine read_data_and_plot_traj


subroutine set_criteria( index, x0_, z0_, px_0_, pz_0_, pinf_x_, pinf_z_, ts_, amp_M_, ierr_, i_type)

    implicit none
    integer:: index
    double precision:: x0_, z0_, px_0_, pz_0_, pinf_x_, pinf_z_
    double complex:: ts_, amp_M_
    integer:: ierr_, i_tyep_
    double precision:: r
    integer:: i_type_, i_type


    if( dreal(ts_) < 560 .and. dreal(ts_) > 260 ) then
        write(*, '(i4, 6(e15.8,1x))'), index, ts_, pinf_x_, pinf_z_, amp_M_;
        print*, ierr_;
        write( FID_FILTER, '(i3)'), index;
        write( FID_TRAJ_M, '(2(e15.8,1x), i4)'), amp_M_, i_type;
    end if

    return;
    
end subroutine set_criteria
