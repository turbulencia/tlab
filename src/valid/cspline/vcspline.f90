#include "types.h"
#include "dns_const.h"

!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2022/03/30 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate cubic splines with different boundary conditions.
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################

program CSPLINE

    use TLAB_TYPES, only: grid_dt

    implicit none

    ! define spline parameters here
    type(grid_dt) :: g, g_int                            ! original and refined grid
    TINTEGER, parameter :: inb_grid = 57
    TINTEGER, parameter :: imax = 11                       ! number of data points
    TINTEGER, parameter :: mesh = 10                       ! mesh refinement factor (mesh=1 for x_int=x)
    TINTEGER, parameter :: imax_int = (imax + (mesh - 1)*(imax - 1)) ! number of spline data points
    TINTEGER :: i
    TREAL :: xb, xe, lambda, r
    TREAL :: res_2, res_inf
    logical :: periodic, random, uniform
    ! data arrays
    TREAL, dimension(imax, inb_grid) :: x
    TREAL, dimension(imax_int, inb_grid) :: x_int
    TREAL, dimension(imax) :: y
    TREAL, dimension(imax_int) :: y_sp, y_int, delta
    TREAL, dimension(imax_int) :: dydx, ddydx
    ! boundary conditions
    TINTEGER, dimension(2) :: bc    ! boundary condition for both endpoints
    TREAL, dimension(2) :: bcval ! boundary values of 1st or 2nd deriv. at endpoints
    ! working arrays
    TREAL, dimension(11*imax) :: wrk
    TREAL, dimension(imax, 5) :: wrk1d
    TREAL, dimension(imax_int, 5) :: wrk1d_int

    integer, parameter :: i0 = 0, i1 = 1

! ###################################################################
! test data (random data or specific test function)
    ! random = .false.  ! uniform grid with sin() function
    random = .true. ! non-uniform grid with random function values

! choose boundary type and boundary values (frist/second derivative)
! for cubic splines here
    ! bc(1) = CS_BCS_CLAMPED;  bcval(:) = C_0_R  ! 1st derv. is zero
    ! bc(1) = CS_BCS_FIXED_1;  bcval(:) = C_3_R  ! fixed 1st deriv.
    ! bc(1) = CS_BCS_NATURAL;  bcval(:) = C_0_R  ! 2nd deriv. is zero
    ! bc(1) = CS_BCS_FIXED_2;  bcval(:) = C_3_R  ! fixed 2nd deriv.
    bc(1) = CS_BCS_PERIODIC; bcval(:) = C_0_R    ! periodic BCs (for start and end of curve!)
    bc(2) = bc(1)                                ! also mixed BCs possible (not for periodic BC!)

! ###################################################################
! initialize original grid with test function
    if (bc(1) == CS_BCS_PERIODIC) then; periodic = .true.; bc(2) = CS_BCS_PERIODIC
    else; periodic = .false.; end if
    if (.not. random) then; uniform = .true.
    else; uniform = .false.; end if

    g%periodic = periodic; g%uniform = uniform
    g%size = imax
    g%scale = C_1_R
    lambda = C_1_R
    g%mode_fdm = FDM_COM6_JACOBIAN

    if (random) then
        do i = 1, imax
            call random_number(r)
            x(i, 1) = i - C_1_R + r
            call random_number(r)
            y(i) = r
        end do
        x(1, 1) = C_0_R
    else
        do i = 1, imax
            x(i, 1) = M_REAL(i - 1)/M_REAL(imax - 1)*g%scale
            y(i) = sin(C_2_R*C_PI_R/g%scale*lambda*x(i, 1))
        end do
    end if
    if (periodic) y(imax) = y(1)

! set interval for spline approximation
    xb = x(1, 1); xe = x(imax, 1)

! create equally spaced array to evaluate spline with boundaries [xb,xe]
    g_int%periodic = g%periodic; g_int%size = imax_int; g_int%scale = g%scale
    g_int%uniform = g%uniform; g_int%mode_fdm = g%mode_fdm
    do i = 1, imax_int
        x_int(i, 1) = xb + (xe - xb)*(i - C_1_R)/(imax_int - C_1_R)
    end do
    if (random) then
        do i = 1, imax_int
            y_int(i) = C_0_R ! exact function exists only for non-random data
        end do
    else
        do i = 1, imax_int
            y_int(i) = sin(C_2_R*C_PI_R/g_int%scale*lambda*x_int(i, 1))
        end do
    end if

! initialize grids for fdm calls
    call FDM_INITIALIZE(x, g, wrk1d)
    call FDM_INITIALIZE(x_int, g_int, wrk1d_int)

! cubic spline function
    call CUBIC_SPLINE(bc, bcval, imax, imax_int, g%nodes, y, g_int%nodes, y_sp, wrk)

! ###################################################################
    write (*, *) '================================================================================'
    write (*, *) '================  Validation routines for cubic splines ========================'
    write (*, *) '================================================================================'
! ###################################################################
! 1. Validation (works only for uniform data)
    write (*, *) '1. Validation routine: Check if data points coincident with spline points'
    if (.not. random) then
        do i = 1, imax
            if ((abs(y(i) - y_sp(1 + (i - 1)*mesh))) <= 1e-10) then
                write (*, 40) i, g%nodes(i), y(i), (1 + (i - 1)*mesh), abs(y(i) - y_sp(1 + (i - 1)*mesh))
            else
        write(*,50) i, g%nodes(i), y(i), (1 + (i-1)*mesh), g_int%nodes(1 + (i-1)*mesh), y_sp(1 + (i-1)*mesh), abs(y(i) - y_sp(1 + (i-1)*mesh))
            end if
        end do
    else
        do i = 1, imax, imax - 1
            if ((abs(y(i) - y_sp(1 + (i - 1)*mesh))) <= 1e-10) then
                write (*, 40) i, g%nodes(i), y(i), (1 + (i - 1)*mesh), abs(y(i) - y_sp(1 + (i - 1)*mesh))
            else
        write(*,50) i, g%nodes(i), y(i), (1 + (i-1)*mesh), g_int%nodes(1 + (i-1)*mesh), y_sp(1 + (i-1)*mesh), abs(y(i) - y_sp(1 + (i-1)*mesh))
            end if
        end do
    end if
40  format(1x, 'Data point ', I3, ' with ', '(', F5.2, ',', F5.2, ')', ' is on spline point ', I3, ', residua = ', ES9.2)
50  format(1X,'Data point ', I3, ' with ', '(', F5.2, ',', F5.2,')', ' is not on spline point ', I3, ' with ', '(', F4.2, ',', F4.2,')',', residua = ' ,ES9.2)

! 2. Validation
! compute first and second derivative at boundary points
    if ((g%mode_fdm == FDM_COM6_JACOBIAN) .and. (.not. periodic)) then
        ! first derivative
        ! call FDM_C1N6_LHS(imax_int, i0, i0, g_int%jac, wrk1d_int(1, 1), wrk1d_int(1, 2), wrk1d_int(1, 3))
        ! call FDM_C1N6_RHS(imax_int, i1, i0, i0, y_sp, dydx)
        call TRIDFS(imax_int, wrk1d_int(1, 1), wrk1d_int(1, 2), wrk1d_int(1, 3))
        call TRIDSS(imax_int, i1, wrk1d_int(1, 1), wrk1d_int(1, 2), wrk1d_int(1, 3), dydx)
        ! second derivative
        ! call FDM_C2N6H_LHS(imax_int, i0, i0, g_int%jac, wrk1d_int(1, 1), wrk1d_int(1, 2), wrk1d_int(1, 3))
        ! call FDM_C2N6H_RHS(imax_int, i1, i0, i0, y_sp, ddydx)
        print*,'to be rewritten using new FDM routines'
        call TRIDFS(imax_int, wrk1d_int(1, 1), wrk1d_int(1, 2), wrk1d_int(1, 3))
        call TRIDSS(imax_int, i1, wrk1d_int(1, 1), wrk1d_int(1, 2), wrk1d_int(1, 3), ddydx)
        write (*, *) '================================================================================'
        write (*, *) '2. Validation routine: derivative values on start end point of curve'
        write (*, '(1X,A,ES9.2,3X,ES9.2)') '1st derivative values at boundary points = ', dydx(1), dydx(imax_int)
        write (*, '(1X,A,ES9.2,3X,ES9.2)') '2nd derivative values at boundary points = ', ddydx(1), ddydx(imax_int)
    end if

! 3. Validation
    if (.not. random) then
        do i = 1, imax_int
            delta(i) = y_int(i) - y_sp(i)
        end do
        res_2 = norm2(delta); res_inf = maxval(abs(delta))
        write (*, *) '================================================================================'
        write (*, *) '3. Validation routine: Error norms of exact solution and spline points'
        write (*, '(1X,A,ES9.2)') 'L_2   error norm                         = ', res_2
        write (*, '(1X,A,ES9.2)') 'L_inf error norm                         = ', res_inf
    end if
    write (*, *) '================================================================================'

! IO - error and function values
    open (20, file='cspline_org.dat')
    do i = 1, imax
        write (20, 1000) g%nodes(i), y(i)
    end do
    close (20)
    !
    open (30, file='cspline_int.dat')
    do i = 1, imax_int
        write (30, 1000) g_int%nodes(i), y_int(i), y_sp(i), y_int(i) - y_sp(i)
    end do
    close (30)
1000 format(6(1x, e16.10))

    ! Validation option in python with:
    !   from   scipy.interpolate import CubicSpline
    !   sp = CubicSpline(xorg,yorg,bc_type=((1,0),(1,0)))

    stop

end program CSPLINE
