subroutine display_info( ts )
    use mod_p0, only: p0_z, p0_x
    implicit none
    double complex, intent(in):: ts

    double complex, external:: sub_traj_x_0, sub_traj_z_0, sub_traj_vx_0, sub_traj_vz_0

    write(*, '(2(a,e15.8,x))'),'p0_z =', p0_z, '  p0_x =', p0_x
    write(*, '(a,2(e15.8,x))'), 'ts =', ts
    write(*, '(a,2(e15.8,x))'), 'sub_x_0(ts):', sub_traj_x_0( ts, ts )
    write(*, '(a,2(e15.8,x))'), 'sub_z_0(ts):', sub_traj_z_0( ts, ts )
    write(*, '(a,2(e15.8,x))'), 'sub_vx_0(ts):', sub_traj_vx_0( ts )
    write(*, '(a,2(e15.8,x))'), 'sub_vz_0(ts):', sub_traj_vz_0( ts )

    return
end subroutine display_info
