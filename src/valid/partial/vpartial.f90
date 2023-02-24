!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/17 - J. Kostelecky
!#              Modified for pentadiagonal schemes
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact schemes.
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

program VPARTIAL

    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: C1N6M_ALPHA
    use OPR_PARTIAL

    implicit none

    type(grid_dt) :: g

    TINTEGER :: jmax, kmax, i, l
    TINTEGER, parameter :: imax = 128, len = 1, inb_grid = 57

    TREAL, dimension(imax, inb_grid) :: x
    TREAL, dimension(len, imax) :: u
    TREAL, dimension(len, imax) :: du1_a, du1_b, du1_c
    TREAL, dimension(len, imax) :: du2_a
    !  TREAL,    DIMENSION(len,imax)      :: du2_n1, du2_n2, du2_n3
    TREAL, dimension(imax, 7) :: wrk1d
    TREAL, dimension(len, 2) :: bcs
    TINTEGER bcs_aux(2, 2)
    TREAL :: lambda, error, dummy
    TINTEGER :: test_type, ibc

    integer, parameter :: i1 = 1

! ###################################################################
! Initialize
    g%size = imax
    g%scale = C_1_R
    g%uniform = .true.
    jmax = 1
    kmax = 1

! Valid stettings
    g%periodic = .false.
    lambda = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda
    test_type = 2
    ibc = 3
    g%mode_fdm = FDM_COM6_JACOBIAN ! FDM_COM6_JACPENTA

    if (g%mode_fdm == FDM_COM6_JACOBIAN) C1N6M_ALPHA = 0.56

!  ###################################################################

    if (g%periodic) then
        do i = 1, imax
            x(i, 1) = M_REAL(i - 1)/M_REAL(imax)*g%scale
        end do
    else
        do i = 1, imax
            x(i, 1) = M_REAL(i - 1)/M_REAL(imax - 1)*g%scale
        end do
    end if

    call FDM_INITIALIZE(x, g, wrk1d)

! Bcs
    bcs(1, 1) = C_0_R
    bcs(1, 2) = C_0_R
    bcs_aux = 0

! ###################################################################
! Define the function and analytic derivatives
    do i = 1, imax
        do l = 1, len
! single-mode
            u(l, i) = &
                sin(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))!+C_PI_R/C_4_R)
            du1_a(l, i) = (C_2_R*C_PI_R/g%scale*lambda) &
                          *cos(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))!+C_PI_R/C_4_R)
            du2_a(l, i) = -(C_2_R*C_PI_R/g%scale*lambda)**2 &
                          *u(l, i)
! Gaussian
            ! dummy = C_1_R / ( C_2_R*(g%scale/M_REAL(lambda*l))**2 )
            ! u(l,i)     = EXP(-dummy*(g%nodes(i)-x_0*g%scale)**2)
            ! du1_a(l,i) =-C_2_R *dummy *(g%nodes(i)-x_0*g%scale) *u(l,i)
            ! du2_a(l,i) =-C_2_R *dummy *(g%nodes(i)-x_0*g%scale) *du1_a(l,i) - C_2_R *dummy *u(l,i)
! Exponential
            ! u(l,i)     = EXP(g%nodes(i)/(g%scale/lambda))
            ! du1_a(l,i) = lambda/g%scale*u(l,i)
            ! du2_a(l,i) = lambda/g%scale*du1_a(l,i)
! delta-function
            ! u(i)     = MAX(C_0_R,C_2_R-M_REAL(i))
            ! du1_a(i) = C_0_R
            ! du2_a(i) = C_0_R
! hyperboic tangent
            ! u(l,i)     = lambda*LOG(C_1_R+EXP(g%nodes(i)/lambda))
            ! du1_a(l,i) = C_05_R*(C_1_R+TANH(C_05_R*g%nodes(i)/lambda))
            ! du2_a(l,i) = C_025_R/lambda/(COSH(C_05_R*g%nodes(i)/lambda))**2
! Polynomial
            ! dummy = C_4_R
            ! u(l,i)     =                       ( (g%scale-g%nodes(i)) /lambda)** dummy
            ! du1_a(l,i) = dummy                *( (g%scale-g%nodes(i)) /lambda)**(dummy-C_1_R)
            ! du2_a(l,i) = dummy *(dummy-C_1_R) *( (g%scale-g%nodes(i)) /lambda)**(dummy-C_2_R)
        end do
    end do

! ###################################################################
! Testing first-order derivatives

    if (test_type == 1) then
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, du1_b)

! ! -------------------------------------------------------------------
! ! Testing second-order derivatives
! ! -------------------------------------------------------------------
!   ! Jacobian based
!     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g, u,     du2_n2, du1_n)
!     CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g, du1_n, du2_n1)
!   ! Direct metrics
!     CALL FDM_C2N6ND_INITIALIZE(imax, x, wrk1d(1,1), wrk1d(1,4))
!     CALL TRIDFS(imax,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))
!     CALL FDM_C2N6ND_RHS(imax,len, wrk1d(1,4), u, du2_n3)
!     CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du2_n3)

! -------------------------------------------------------------------
    elseif (test_type == 2) then ! Testing new BCs routines

        if (g%mode_fdm == FDM_COM6_JACOBIAN) then
            call FDM_C1N6_BCS_LHS(imax, ibc, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1N6_BCS_RHS(imax, len, ibc, u, du1_b)
        elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
            call FDM_C1N6M_BCS_LHS(imax, ibc, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call FDM_C1N6M_BCS_RHS(imax, len, ibc, u, du1_b)
        end if

        if (ibc == 0) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
                call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du1_b)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                call PENTADFS2(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
                call PENTADSS2(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), du1_b)
            end if

        else if (ibc == 1) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                wrk1d(:, 4) = C_0_R
                wrk1d(1, 4) = C_1_R
                wrk1d(2, 4) = wrk1d(1, 1)
                wrk1d(3, 4) = wrk1d(2, 1)
                call TRIDFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3))
                call TRIDSS(imax - 1, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), du1_b(1, 2))
                call TRIDSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 4) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 2)*du1_b(:, 1) + wrk1d(1, 3)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                wrk1d(:, 6) = C_0_R
                wrk1d(1, 6) = C_1_R
                wrk1d(2, 6) = wrk1d(1, 2)
                wrk1d(3, 6) = wrk1d(2, 2)
                call PENTADFS2(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
                call PENTADSS2(imax - 1, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), du1_b(1, 2))
                call PENTADSS2(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 6) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 3)*du1_b(:, 1) + wrk1d(1, 4)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
            end if

        else if (ibc == 2) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                wrk1d(:, 5) = C_0_R
                wrk1d(imax, 5) = C_1_R
                wrk1d(imax - 2, 5) = wrk1d(imax - 1, 3)
                wrk1d(imax - 1, 5) = wrk1d(imax, 3)
                call TRIDFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
                call TRIDSS(imax - 1, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du1_b)
                call TRIDSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 5))
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 5) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 1)*du1_b(:, imax - 1) + wrk1d(imax, 2)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                wrk1d(:, 6) = C_0_R
                wrk1d(imax, 6) = C_1_R
                wrk1d(imax - 2, 6) = wrk1d(imax - 1, 4)
                wrk1d(imax - 1, 6) = wrk1d(imax, 4)
                call PENTADFS2(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
                call PENTADSS2(imax - 1, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), du1_b)
                call PENTADSS2(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 6) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 2)*du1_b(:, imax - 1) + wrk1d(imax, 3)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            end if

        else if (ibc == 3) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                wrk1d(:, 4) = C_0_R
                wrk1d(1, 4) = C_1_R
                wrk1d(2, 4) = wrk1d(1, 1)
                wrk1d(3, 4) = wrk1d(2, 1)
                wrk1d(:, 5) = C_0_R
                wrk1d(imax, 5) = C_1_R
                wrk1d(imax - 2, 5) = wrk1d(imax - 1, 3)
                wrk1d(imax - 1, 5) = wrk1d(imax, 3)
                call TRIDFS(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3))
                call TRIDSS(imax - 2, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), du1_b(1, 2))
                call TRIDSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4))
                call TRIDSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 5))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 4) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 2)*du1_b(:, 1) + wrk1d(1, 3)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 5) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 1)*du1_b(:, imax - 1) + wrk1d(imax, 2)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                wrk1d(:, 6) = C_0_R
                wrk1d(1, 6) = C_1_R
                wrk1d(2, 6) = wrk1d(1, 2)
                wrk1d(3, 6) = wrk1d(2, 2)
                wrk1d(:, 7) = C_0_R
                wrk1d(imax, 7) = C_1_R
                wrk1d(imax - 2, 7) = wrk1d(imax - 1, 4)
                wrk1d(imax - 1, 7) = wrk1d(imax, 4)
                call PENTADFS2(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
                call PENTADSS2(imax - 2, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), du1_b(1, 2))
                call PENTADSS2(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                call PENTADSS2(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 6) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 3)*du1_b(:, 1) + wrk1d(1, 4)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = C_0_R
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 7) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 2)*du1_b(:, imax - 1) + wrk1d(imax, 3)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            end if
        end if
    end if

! ###################################################################
! IO - Error and function values
    open (20, file='partial.dat')
    error = C_0_R; dummy = C_0_R
    do i = 1, imax
        do l = 1, len
            ! Testing first-order derivatives
            write (20, 1000) g%nodes(i), u(l, i), du1_a(l, i), du1_b(l, i), du1_a(l, i) - du1_b(l, i)
            du1_c(l, i) = abs(du1_a(l, i) - du1_b(l, i))
            dummy = dummy + du1_a(l, i)*du1_a(l, i)
            error = error + du1_c(l, i)*du1_c(l, i)
            ! ! Testing second-order derivatives
            ! WRITE(20,1000) g%nodes(i), u(l,i), du2_a(l,i), du2_n2(l,i), du2_a(l,i)-du2_n2(l,i)
            ! du1_c(l,i)= ABS(du2_a(l,i)-du2_n2(l,i))
            ! dummy = dummy + du2_a(l,i)*du2_a(l,i)
            ! error = error + du1_c(l,i) * du1_c(l,i)
        end do
    end do
    close (20)

    write (*, *) 'Solution L2-norm ...........:', sqrt(g%jac(1, 1)*dummy)/M_REAL(len)
    if (dummy == C_0_R) stop
    write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error)/maxval(abs(du1_a))
    write (*, *) 'Relative Error Linf-norm ...:', maxval(du1_c(1, 1:imax))/maxval(abs(du1_a))

    stop

1000 format(5(1x, e12.5))

end program VPARTIAL