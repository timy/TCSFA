#include '../../include/inc_rk4.h'

subroutine rk4_sub( h, t0, y0, y, newton_equation )

    implicit none
    integer, parameter:: ne = RK4_NE
    double precision:: h, t0, t
    double complex:: y0(ne), y(ne)
    double complex:: k1(ne), k2(ne), k3(ne), k4(ne)
    external:: newton_equation

    t = t0
    y = y0

    call newton_equation( ne, t, y, k1 )

    t = t0 + 0.5d0 * h

    y = y0 + 0.5d0 * h * k1
    call newton_equation( ne, t, y, k2 )

    y = y0 + 0.5d0 * h * k2
    call newton_equation( ne, t, y, k3 )

    t = t0 + h
    y = y0 + h * k3
    call newton_equation( ne, t, y, k4 )

    y = y0 + h * ( k1 + 2d0 * ( k2 + k3 ) + k4 ) / 6d0

    return
end subroutine rk4_sub
