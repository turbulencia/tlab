!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/20 - J. Kostelecky
!#              Modified for pentadiagonal schemes
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate.
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

program VINTEGRAL

    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: C1N6M_ALPHA
    use OPR_PARTIAL

    implicit none

    type(grid_dt) :: g
    TINTEGER :: jmax, kmax, i
    TINTEGER, parameter :: imax = 128, inb_grid = 57
    TREAL, dimension(imax, inb_grid) :: x
    TREAL, dimension(imax) :: u, w_n, f, tmp
    TREAL, dimension(imax) :: du1_a, dw1_n, du2_a
    TREAL, dimension(imax, 8) :: wrk1d
    integer, dimension(2, 2) :: bcs
    TREAL :: lambda, error, sol, dummy, wk, x_0
    TINTEGER :: test_type, ibc
    TINTEGER :: imin_loc, imax_loc

    integer, parameter :: i1 = 1
    
! ###################################################################
! Initialize
    g%size = imax
    g%scale = C_1_R
    g%uniform = .true.
    jmax = 1
    kmax = 1
    bcs = 0

! Valid stettings
    g%periodic = .true.
    wk = 1 ! WRITE(*,*) 'Wavenumber ?'; READ(*,*) wk
    lambda = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda
    test_type = 1
    g%mode_fdm = FDM_COM6_JACOBIAN ! FDM_COM6_JACPENTA

    if (g%mode_fdm == FDM_COM6_JACPENTA) C1N6M_ALPHA = 0.56

! ###################################################################

    do i = 1, imax
        x(i, 1) = M_REAL(i - 1)/M_REAL(imax - 1)*g%scale
    end do

    call FDM_INITIALIZE(x, g, wrk1d)

    x_0 = C_075_R

    wrk1d = C_0_R

! ###################################################################
! Define the function f and analytic derivatives
    do i = 1, imax
! single-mode
        u(i) = sin(C_2_R*C_PI_R/g%scale*wk*g%nodes(i) + C_PI_R/C_4_R)
        du1_a(i) = (C_2_R*C_PI_R/g%scale*wk) &
                   *cos(C_2_R*C_PI_R/g%scale*wk*g%nodes(i) + C_PI_R/C_4_R)
        ! u(i)     =              COS(C_2_R*C_PI_R*g%nodes(i))
        ! du1_a(i) =-C_PI_R*C_2_R*SIN(C_2_R*C_PI_R*g%nodes(i))
! Gaussian
        ! u(i)     = EXP(-(g%nodes(i)-x_0*g%scale)**2/(C_2_R*(g%scale/wk)**2))
        ! du1_a(i) =-(g%nodes(i)-x_0*g%scale)/(g%scale/wk)**2*u(i)
        ! v(i)     = EXP(-(g%nodes(i)-x_0*C_05_R*g%scale)**2/(C_2_R*(g%scale/wk)**2))
        ! dv1_a(i) =-(g%nodes(i)-x_0*C_05_R*g%scale)/(g%scale/wk)**2*v(i)
! exponential
        ! u(i) = EXP(-g%nodes(i)*lambda)
        ! du1_a(i) = -lambda*u(i)
        ! du2_a(i) =  lambda*lambda*u(i)
        ! v(i) = EXP(g%nodes(i)*x_0/g%scale)
        ! dv1_a(i) = x_0/g%scale*v(i)
! step
        ! u(i) = MAX(C_0_R,(g%nodes(i)-g%nodes(imax/2))*x_0)
        ! du1_a(i) = (C_1_R+SIGN(C_1_R,g%nodes(i)-g%nodes(imax/2)))*C_05_R*x_0
! tanh
        ! u(i) = x_0 * LOG( C_1_R + EXP( (g%nodes(i)-g%nodes(imax/2))/x_0 ) )
        ! du1_a(i) = C_05_R*( C_1_R + TANH( C_05_R*(g%nodes(i)-g%nodes(imax/2))/x_0 ) )
! polynomial
        ! u(i)     =         g%nodes(i)** lambda
        ! du1_a(i) = lambda*(g%nodes(i)**(lambda-C_1_R))
! zero
        ! u(i) = C_0_R
        ! du1_a(i) = C_0_R
    end do

! ###################################################################
! Integral
! ###################################################################
    if (test_type == 0) then

        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g, u, f) ! f = du1_a

        ibc = 2

        if (g%mode_fdm == FDM_COM6_JACOBIAN) then
            call INT_C1N6_LHS(imax, ibc, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call INT_C1N6_RHS(imax, i1, ibc, g%jac, f, w_n)
        elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
            call INT_C1N6M_LHS(imax, ibc, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
            call INT_C1N6M_RHS(imax, i1, ibc, g%jac, f, w_n)
        end if

        if (ibc == 1) then ! at the bottom
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                call PENTADFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
                call PENTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), w_n(2))
                w_n(1) = C_0_R; w_n = w_n + u(1)    ! BCs
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                call HEPTADFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7))
                call HEPTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7), w_n(2))
                w_n(1) = C_0_R; w_n = w_n + u(1)    ! BCs
            end if

        else if (ibc == 2) then ! at the top
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                call PENTADFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
                call PENTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), w_n(1))
                w_n(imax) = C_0_R; w_n = w_n + u(imax) ! BCs
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                call HEPTADFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                call HEPTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7), w_n(1))
                w_n(imax) = C_0_R; w_n = w_n + u(imax) ! BCs
            end if
        end if

! ###################################################################
! First order equation
! ###################################################################
    else if (test_type == 1) then

        f = du1_a + lambda*u

        ibc = 2

        if (g%mode_fdm == FDM_COM6_JACOBIAN) then
            call INT_C1N6_LHS_E(imax, ibc, g%jac, lambda, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
            call INT_C1N6_RHS(imax, i1, ibc, g%jac, f, w_n)
        elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
            call INT_C1N6M_LHS_E(imax,    ibc, g%jac, lambda, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),wrk1d(1,6),wrk1d(1,7), wrk1d(1,8))
            call INT_C1N6M_RHS(imax, i1, ibc, g%jac, f, w_n)
        end if

        if (ibc == 1) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                call PENTADFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
                call PENTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), w_n(2))
                call PENTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                dummy = w_n(1); w_n(1) = C_0_R
                w_n = w_n + u(1)*wrk1d(1:imax, 6) ! BCs
                dummy = (dummy + wrk1d(1, 3)*w_n(1) + wrk1d(1, 4)*w_n(2) + wrk1d(1, 5)*w_n(3))/g%jac(1, 1)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                call HEPTADFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7))
                call HEPTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7), w_n(2))
                call HEPTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7), wrk1d(2, 8))
                dummy = w_n(1); w_n(1) = C_0_R
                w_n = w_n + u(1)*wrk1d(1:imax, 8) ! BCs
                dummy = (dummy + wrk1d(1, 4)*w_n(1) + wrk1d(1, 5)*w_n(2) + wrk1d(1, 6)*w_n(3))/g%jac(1, 1)
            end if

        else if (ibc == 2) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                call PENTADFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
                call PENTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), w_n(1))
                call PENTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
                dummy = w_n(imax); w_n(imax) = C_0_R
                w_n = w_n + u(imax)*wrk1d(1:imax, 6) ! BCs
             dummy = (dummy + wrk1d(imax, 1)*w_n(imax - 2) + wrk1d(imax, 2)*w_n(imax - 1) + wrk1d(imax, 3)*w_n(imax))/g%jac(imax, 1)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                call HEPTADFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                call HEPTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7), w_n(1))
                call HEPTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7), wrk1d(1, 8))
                dummy = w_n(imax); w_n(imax) = C_0_R
                w_n = w_n + u(imax)*wrk1d(1:imax, 8) ! BCs
             dummy = (dummy + wrk1d(imax, 2)*w_n(imax - 2) + wrk1d(imax, 3)*w_n(imax - 1) + wrk1d(imax, 4)*w_n(imax))/g%jac(imax, 1)
            end if
        end if

        if (ibc == 1) then
            write (*, *) dummy, dw1_n(1)
        else
            write (*, *) dummy, dw1_n(imax)
        end if

        f = u

! ###################################################################
! Second order equation
! ###################################################################
    else if (test_type == 2) then

        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g, u, f, tmp)

        dummy = lambda*lambda
        f = du2_a - dummy*u

        call INT_C2N6_LHS_E(imax, g%jac, dummy, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
        call PENTADFS(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
        call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
        call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
        call INT_C2N6_RHS(imax, i1, g%jac, f, w_n)
        call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), w_n(2))
        w_n(:) = w_n(:) + u(1)*wrk1d(:, 6) + u(imax)*wrk1d(:, 7)

    end if

! ###################################################################
! IO - Error and function values
    open (20, file='integral.dat')
    error = C_0_R
    sol = C_0_R
    imin_loc = 1; imax_loc = imax
    do i = imin_loc, imax_loc
        write (20, 1000) g%nodes(i), u(i), w_n(i), u(i) - w_n(i)
        w_n(i) = abs(u(i) - w_n(i))
        error = error + w_n(i)*w_n(i)
        sol = sol + u(i)*u(i)
    end do
    close (20)

    write (*, *) 'Solution L2-norm ......:', sqrt(g%jac(1, 1)*sol)
    write (*, *) 'Error L2-norm .........:', sqrt(g%jac(1, 1)*error)
    write (*, *) 'Error Linf-norm .......:', maxval(w_n(1:imax))
    write (*, *) 'Relative error ........:', sqrt(error)/sqrt(sol)
    write (*, *) 'Derivative overshoot ..:', minval(dw1_n(1:imax))

    stop

1000 format(6(1x, e17.10e3))

end program VINTEGRAL
