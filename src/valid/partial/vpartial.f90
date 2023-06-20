#include "dns_const.h"

program VPARTIAL
    use TLAB_CONSTANTS
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk3d, inb_txc, isize_txc_field
    use TLAB_VARS, only: reynolds, schmidt
    ! use TLAB_VARS, only: C1N6M_ALPHA
    use TLAB_PROCS
    use TLAB_ARRAYS, only: wrk1d, txc, x
    use FDM_COM_DIRECT
    use FDM_PROCS, only: MatMul_5d
    use OPR_PARTIAL

    implicit none

    type(grid_dt) :: g

    integer(wi) :: i, l, len

    real(wp), dimension(:, :), pointer :: u
    real(wp), dimension(:, :), pointer :: du1_a, du1_b, du1_c, du1_n
    real(wp), dimension(:, :), pointer :: du2_a, du2_n1, du2_n2, du2_n3
    real(wp), dimension(:, :), pointer :: bcs
    integer(wi) bcs_aux(2, 2)
    real(wp) :: lambda, error, dummy
    integer(wi) :: test_type, ibc

    integer, parameter :: i1 = 1

! ###################################################################
! Initialize
    imax = 129
    jmax = 1
    kmax = 1
    len = jmax*kmax

    reynolds = 1.0_wp   ! Needed in FDM_INITIALIZE
    schmidt = 1.0_wp

    g%inb_grid = 57
    g%size = imax
    g%scale = 1.0_wp
    g%uniform = .false.

    isize_field = imax*jmax*kmax
    isize_txc_field = isize_field
    isize_wrk3d = isize_txc_field
    isize_wrk1d = imax
    inb_wrk1d = 20
    inb_txc = 9

    call TLAB_ALLOCATE(__FILE__)

    u(1:len, 1:imax) => txc(1:imax*jmax*kmax, 1)
    du1_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 2)
    du1_b(1:len, 1:imax) => txc(1:imax*jmax*kmax, 3)
    du1_c(1:len, 1:imax) => txc(1:imax*jmax*kmax, 4)
    du1_n(1:len, 1:imax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 6)
    du2_n1(1:len, 1:imax) => txc(1:imax*jmax*kmax, 7)
    du2_n2(1:len, 1:imax) => txc(1:imax*jmax*kmax, 8)
    du2_n3(1:len, 1:imax) => txc(1:imax*jmax*kmax, 9)

    allocate (bcs(len, 2))
    call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, x, [g%size, g%inb_grid], g%name)

    ! Valid settings
    test_type = 1

    g%periodic = .false.
    lambda = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda
    ibc = 3
    g%mode_fdm = FDM_COM6_JACOBIAN ! FDM_COM6_JACPENTA
    ! g%mode_fdm = FDM_COM6_DIRECT

    ! if (g%mode_fdm == FDM_COM6_JACOBIAN) C1N6M_ALPHA = 0.56_wp

!  ###################################################################

    if (g%periodic) then
        do i = 1, imax
            x(i, 1) = real(i - 1, wp)/real(imax, wp)*g%scale
        end do
    else
        ! do i = 1, imax
        !     x(i, 1) = real(i - 1, wp)/real(imax - 1, wp)*g%scale
        ! end do
        open (21, file='y.dat')
        do i = 1, imax
        ! do i = imax, 1, -1
            read (21, *) x(i, 1)
        end do
        close (21)
        ! x(:,1) = x(1,1) -x(:,1)
        g%scale = x(imax, 1) - x(1, 1)
    end if

    call FDM_INITIALIZE(x, g, wrk1d)

! Bcs
    bcs(1, 1) = 0.0_wp
    bcs(1, 2) = 0.0_wp
    bcs_aux = 0

! ###################################################################
! Define the function and analytic derivatives
    do i = 1, imax
        do l = 1, len
! single-mode
            u(l, i) = &
                sin(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))!+pi_wp/C_4_R)
            du1_a(l, i) = (2.0_wp*pi_wp/g%scale*lambda) &
                          *cos(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))!+pi_wp/C_4_R)
            du2_a(l, i) = -(2.0_wp*pi_wp/g%scale*lambda)**2 &
                          *u(l, i)
! Gaussian
            ! dummy = 1.0_wp / ( 2.0_wp*(g%scale/M_REAL(lambda*l))**2 )
            ! u(l,i)     = EXP(-dummy*(g%nodes(i)-x_0*g%scale)**2)
            ! du1_a(l,i) =-2.0_wp *dummy *(g%nodes(i)-x_0*g%scale) *u(l,i)
            ! du2_a(l,i) =-2.0_wp *dummy *(g%nodes(i)-x_0*g%scale) *du1_a(l,i) - 2.0_wp *dummy *u(l,i)
! Exponential
            ! u(l,i)     = EXP(g%nodes(i)/(g%scale/lambda))
            ! du1_a(l,i) = lambda/g%scale*u(l,i)
            ! du2_a(l,i) = lambda/g%scale*du1_a(l,i)
! delta-function
            ! u(i)     = MAX(0.0_wp,2.0_wp-M_REAL(i))
            ! du1_a(i) = 0.0_wp
            ! du2_a(i) = 0.0_wp
! hyperboic tangent
            ! u(l,i)     = lambda*LOG(1.0_wp+EXP(g%nodes(i)/lambda))
            ! du1_a(l,i) = C_05_R*(1.0_wp+TANH(C_05_R*g%nodes(i)/lambda))
            ! du2_a(l,i) = C_025_R/lambda/(COSH(C_05_R*g%nodes(i)/lambda))**2
! Polynomial
            ! dummy = 4.0_wp
            ! u(l,i)     =                       ( (g%scale-g%nodes(i)) /lambda)** dummy
            ! du1_a(l,i) = dummy                *( (g%scale-g%nodes(i)) /lambda)**(dummy-1.0_wp)
            ! du2_a(l,i) = dummy *(dummy-1.0_wp) *( (g%scale-g%nodes(i)) /lambda)**(dummy-2.0_wp)
        end do
    end do

! ###################################################################

    if (test_type == 1) then
        bcs_aux = 0
! Testing first-order derivatives
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, du1_b)

! Testing second-order derivatives
        ! Jacobian based
        ! call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_aux, g, u, du2_n2, du1_n)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, du1_n, du2_n1)
        ! Direct metrics
        CALL FDM_C2N6ND_INITIALIZE(imax, x, wrk1d(1,1), wrk1d(1,4))
        ! CALL FDM_C2N4ND_INITIALIZE(imax, x, wrk1d(1,1), wrk1d(1,4))
        CALL TRIDFS(imax,     wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))
        CALL MatMul_5d(imax,len, wrk1d(1,4), u, du2_n2)
        CALL TRIDSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), du2_n2)

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
                wrk1d(:, 4) = 0.0_wp
                wrk1d(1, 4) = 1.0_wp
                wrk1d(2, 4) = wrk1d(1, 1)
                wrk1d(3, 4) = wrk1d(2, 1)
                call TRIDFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3))
                call TRIDSS(imax - 1, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), du1_b(1, 2))
                call TRIDSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 4) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 2)*du1_b(:, 1) + wrk1d(1, 3)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                wrk1d(:, 6) = 0.0_wp
                wrk1d(1, 6) = 1.0_wp
                wrk1d(2, 6) = wrk1d(1, 2)
                wrk1d(3, 6) = wrk1d(2, 2)
                call PENTADFS2(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
                call PENTADSS2(imax - 1, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), du1_b(1, 2))
                call PENTADSS2(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 6) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 3)*du1_b(:, 1) + wrk1d(1, 4)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
            end if

        else if (ibc == 2) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                wrk1d(:, 5) = 0.0_wp
                wrk1d(imax, 5) = 1.0_wp
                wrk1d(imax - 2, 5) = wrk1d(imax - 1, 3)
                wrk1d(imax - 1, 5) = wrk1d(imax, 3)
                call TRIDFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
                call TRIDSS(imax - 1, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du1_b)
                call TRIDSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 5))
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 5) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 1)*du1_b(:, imax - 1) + wrk1d(imax, 2)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                wrk1d(:, 6) = 0.0_wp
                wrk1d(imax, 6) = 1.0_wp
                wrk1d(imax - 2, 6) = wrk1d(imax - 1, 4)
                wrk1d(imax - 1, 6) = wrk1d(imax, 4)
                call PENTADFS2(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
                call PENTADSS2(imax - 1, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), du1_b)
                call PENTADSS2(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 6) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 2)*du1_b(:, imax - 1) + wrk1d(imax, 3)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            end if

        else if (ibc == 3) then
            if (g%mode_fdm == FDM_COM6_JACOBIAN) then
                wrk1d(:, 4) = 0.0_wp
                wrk1d(1, 4) = 1.0_wp
                wrk1d(2, 4) = wrk1d(1, 1)
                wrk1d(3, 4) = wrk1d(2, 1)
                wrk1d(:, 5) = 0.0_wp
                wrk1d(imax, 5) = 1.0_wp
                wrk1d(imax - 2, 5) = wrk1d(imax - 1, 3)
                wrk1d(imax - 1, 5) = wrk1d(imax, 3)
                call TRIDFS(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3))
                call TRIDSS(imax - 2, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), du1_b(1, 2))
                call TRIDSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4))
                call TRIDSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 5))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 4) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 2)*du1_b(:, 1) + wrk1d(1, 3)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 5) ! BCs
                end do
                bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 1)*du1_b(:, imax - 1) + wrk1d(imax, 2)*du1_b(:, imax))
                write (*, *) bcs(:, 2), u(:, imax)
            elseif (g%mode_fdm == FDM_COM6_JACPENTA) then
                wrk1d(:, 6) = 0.0_wp
                wrk1d(1, 6) = 1.0_wp
                wrk1d(2, 6) = wrk1d(1, 2)
                wrk1d(3, 6) = wrk1d(2, 2)
                wrk1d(:, 7) = 0.0_wp
                wrk1d(imax, 7) = 1.0_wp
                wrk1d(imax - 2, 7) = wrk1d(imax - 1, 4)
                wrk1d(imax - 1, 7) = wrk1d(imax, 4)
                call PENTADFS2(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
                call PENTADSS2(imax - 2, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), du1_b(1, 2))
                call PENTADSS2(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                call PENTADSS2(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
                bcs(:, 1) = du1_b(:, 1)
                du1_b(:, 1) = 0.0_wp
                do i = 1, imax
                    du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 6) ! BCs
                end do
                bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 3)*du1_b(:, 1) + wrk1d(1, 4)*du1_b(:, 2))
                write (*, *) bcs(:, 1), u(:, 1)
                bcs(:, 2) = du1_b(:, imax)
                du1_b(:, imax) = 0.0_wp
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
    error = 0.0_wp; dummy = 0.0_wp
    do i = 1, imax
        do l = 1, len
            ! Testing first-order derivatives
            ! write (20, 1000) g%nodes(i), u(l, i), du1_a(l, i), du1_b(l, i), du1_a(l, i) - du1_b(l, i)
            ! du1_c(l, i) = abs(du1_a(l, i) - du1_b(l, i))
            ! dummy = dummy + du1_a(l, i)*du1_a(l, i)
            ! error = error + du1_c(l, i)*du1_c(l, i)
            ! Testing second-order derivatives
            write (20, 1000) g%nodes(i), u(l, i), du2_a(l, i), du2_n2(l, i), du2_a(l, i) - du2_n2(l, i)
            du1_c(l, i) = abs(du2_a(l, i) - du2_n2(l, i))
            dummy = dummy + du2_a(l, i)*du2_a(l, i)
            error = error + du1_c(l, i)*du1_c(l, i)
        end do
    end do
    close (20)

    write (*, *) 'Solution L2-norm ...........:', sqrt(g%jac(1, 1)*dummy)/real(len, wp)
    if (dummy == 0.0_wp) stop
    write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error)/maxval(abs(du1_a))
    write (*, *) 'Relative Error Linf-norm ...:', maxval(du1_c(1, 1:imax))/maxval(abs(du1_a))

    stop

1000 format(5(1x, e12.5))

end program VPARTIAL
