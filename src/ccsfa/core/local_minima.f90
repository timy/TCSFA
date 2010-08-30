
#include '../include/inc_ts_guess.h'

#define LMS_DEBUG
!#define LMS_PLOT
!#define LMS_FILE_ID 101
!#define LMS_FILE_NAME1 'dat/lms_map.dat'

subroutine local_minima( xa, xb, ya, yb, func, count, guess )

    implicit none;
    integer, parameter:: nx = LMS_NX
    integer, parameter:: ny = LMS_NY
    double precision, intent(in):: xa, xb, ya, yb    
    double precision:: x(nx), y(ny), val(nx,ny)
    double precision:: dx, dy
    integer:: i, j
    double complex:: grid(nx,ny), cw
    double complex, external:: func;
    integer:: count;
    double complex:: guess(LMS_MAX_COUNT)

#ifdef LMS_DEBUG
    if( xa > xb ) stop 'local_minima error: xa > xb';
    if( ya > yb ) stop 'local_minima error: ya > yb'
#endif

    dx = ( xb - xa ) / ( nx - 1d0 );
    dy = ( yb - ya ) / ( ny - 1d0 );
    forall(i=1:nx) x(i) = xa + dx * ( i - 1d0 );
    forall(j=1:ny) y(j) = ya + dy * ( j - 1d0 );
    
    forall(i=1:nx,j=1:ny) grid(i,j) = dcmplx( x(i), y(j) );


    do j = 1, ny
        do i = 1, nx
            val(i,j) = cdabs( func( dcmplx( x(i), y(j) ) ) );
        end do
    end do

!!$#ifdef LMS_PLOT
!!$    open( LMS_FILE_ID, file=LMS_FILE_NAME1 )
!!$    do j = 1, ny
!!$        do i = 1, nx
!!$            write(LMS_FILE_ID, '(3(e15.8,2x))'), x(i), y(j), val(i,j);
!!$        end do
!!$    end do
!!$    close( LMS_FILE_ID );
!!$#endif

    count = 0;
    do j = 2, ny-1
        do i = 2, nx-1
            if (  val(i,j) <= val(i-1,j) .and. &
                  val(i,j) <= val(i+1,j) .and. &
                  val(i,j) <= val(i,j-1) .and. &
                  val(i,j) <= val(i,j+1) ) then
                count = count + 1;
                if(count > LMS_MAX_COUNT) stop 'exceed LMS_MAX_COUNT';
                guess(count) = dcmplx( x(i), y(j) );
            end if
        end do
    end do

    return;
end subroutine local_minima
