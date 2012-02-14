double complex function trapezoid_sub(a, b, n_step, ts, f)
    implicit none
    double precision, intent(in):: a, b
    integer, intent(in):: n_step
    double complex, intent(in):: ts
    double complex, external:: f
    double precision:: dt
    integer:: i
    
    trapezoid_sub = (0d0, 0d0)
    dt = ( b - a ) / n_step;
    do i = 1, n_step
        trapezoid_sub = trapezoid_sub + f(a+(i-0.5)*dt, ts) * dt
    end do

    return;
end function trapezoid_sub

! simpson rule for sub-barrier correction
double complex function simpson_sub(a, b, n_step, ts, f)
    implicit none
    double precision, intent(in):: a, b
    integer, intent(in):: n_step
    double complex, intent(in):: ts
    double complex, external:: f
    double precision:: dt, h
    integer:: i, m
    m = n_step / 2
    h = (b-a) / (2*m)
    simpson_sub = f(a,ts) + f(b,ts)
    do i = 1, m-1
        simpson_sub = simpson_sub+2d0*f(a+2*i*h,ts)
    end do
    do i = 1, m
        simpson_sub = simpson_sub+4d0*f(a+(2*i-1)*h,ts)
    end do
    simpson_sub = simpson_sub * h / 3d0
    
    return
end function simpson_sub
