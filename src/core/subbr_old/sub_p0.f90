module mod_p0
    implicit none;
    double precision:: p0_z;
    double precision:: p0_x;
end module mod_p0

subroutine set_p0( p0_x_, p0_z_ )
    use mod_p0
    implicit none
    double precision, intent(in):: p0_x_, p0_z_

    p0_x = p0_x_;
    p0_z = p0_z_;
    
    return;
end subroutine set_p0
