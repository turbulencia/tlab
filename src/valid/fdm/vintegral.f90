#include "dns_const.h"

program VINTEGRAL
    use TLAB_CONSTANTS
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLAB_VARS, only: visc, schmidt
    use TLAB_PROCS
    use TLAB_ARRAYS, only: wrk1d, txc, x!, wrk2d, wrk3d
    use FDM_ComX_Direct
    use FDM_Integrate
    use FDM_PROCS
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use OPR_PARTIAL

    implicit none

    type(grid_dt) :: g

    integer(wi) :: i, l, len

    real(wp), dimension(:, :), pointer :: u, w_n, f
    real(wp), dimension(:, :), pointer :: du1_a, dw1_n, du2_a
    real(wp), dimension(:, :), pointer :: bcs
    integer(wi) bcs_aux(2, 2), idr
    real(wp) :: lambda, coef(5), dummy, wk, x_0
    integer(wi) :: test_type, ibc, ip
    integer(wi) :: nmin, nmax, nsize

    integer, parameter :: i1 = 1, cases(2) = [1, 2], cases_new(2) = [BCS_MIN, BCS_MAX]
    real(wp), dimension(:, :), allocatable :: lhs_int, rhs_int

! ###################################################################
! Initialize
    imax = 129
    jmax = 1
    kmax = 1
    len = jmax*kmax

    visc = 1.0_wp   ! Needed in FDM_INITIALIZE
    schmidt = 1.0_wp

    g%inb_grid = 71
    g%size = imax
    g%scale = 1.0_wp
    g%uniform = .false.

    isize_field = imax*jmax*kmax
    isize_txc_field = isize_field
    isize_wrk3d = isize_txc_field
    isize_wrk1d = imax
    isize_wrk2d = len
    inb_wrk1d = 20
    inb_wrk2d = 2
    inb_txc = 9

    call TLAB_ALLOCATE(__FILE__)

    u(1:len, 1:imax) => txc(1:imax*jmax*kmax, 1)
    w_n(1:len, 1:imax) => txc(1:imax*jmax*kmax, 2)
    f(1:len, 1:imax) => txc(1:imax*jmax*kmax, 3)

    du1_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 4)
    dw1_n(1:len, 1:imax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 6)

    allocate (bcs(len, 2))
    call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, x, [g%size, g%inb_grid], g%name)

    allocate (lhs_int(imax, 9), rhs_int(imax, 7))

    g%periodic = .false.
    g%mode_fdm1 = FDM_COM6_JACOBIAN ! FDM_COM6_JACOBIAN_PENTA
    if (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) C1N6M_ALPHA = 0.56
    g%mode_fdm2 = g%mode_fdm1

    wk = 1.0_wp ! WRITE(*,*) 'Wavenumber ?'; READ(*,*) wk
    write (*, *) 'Eigenvalue ?'
    read (*, *) lambda

    test_type = 1

    ! ###################################################################

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
            read (21, *) x(i, 1)
        end do
        close (21)
        g%scale = x(imax, 1) - x(1, 1)
    end if

    call FDM_INITIALIZE(x, g, wrk1d)

    bcs_aux = 0

! ###################################################################
! Define the function f and analytic derivatives
    x_0 = 0.75_wp

    do i = 1, imax
! single-mode
        ! u(:, i) = 1.0_wp + sin(2.0_wp*pi_wp/g%scale*wk*g%nodes(i)) ! + pi_wp/4.0_wp)
        ! du1_a(:, i) = (2.0_wp*pi_wp/g%scale*wk) &
        !               *cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))! + pi_wp/4.0_wp)
! Gaussian
        u(:, i) = exp(-(g%nodes(i) - x_0*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        du1_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*u(:, i)
! exponential
        ! u(:, i) = exp(-g%nodes(i)*wk)
        ! du1_a(:, i) = -wk*u(:, i)
! step
        ! u(:, i) = max(0.0_wp, (g%nodes(i) - g%nodes(imax/2))*x_0)
        ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, g%nodes(i) - g%nodes(imax/2)))*0.5_wp*x_0
! tanh
        ! u(:, i) = x_0*log(1.0_wp + exp((g%nodes(i) - g%nodes(imax/2))/x_0))
        ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(g%nodes(i) - g%nodes(imax/2))/x_0))
! Polynomial
        ! dummy = 4.0_wp
        ! u(:, i) = ((g%scale - g%nodes(i))/wk)**dummy
        ! du1_a(:, i) = -dummy*((g%scale - g%nodes(i))/wk)**(dummy - 1.0_wp)
! zero
        ! u(i) = 0.0_wp
        ! du1_a(i) = 0.0_wp
    end do

! ###################################################################
! First order equation
! ###################################################################
    select case (test_type)
    case (1)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, f)
        f = du1_a + lambda*u

        ! print *, '1. order, Jacobian 4'
        ! call FDM_C1N4_Jacobian(imax, g%jac, g%lu1(:, :), g%rhs1(:, :), coef, g%periodic)
        ! g%nb_diag_1 = [3, 3]
        ! g%rhs1(:, 5) = 0.0_wp
        ! g%rhs1(:, 4) = g%rhs1(:, 3)
        ! g%rhs1(:, 3) = g%rhs1(:, 2)
        ! g%rhs1(:, 2) = g%rhs1(:, 1)
        ! g%rhs1(:, 1) = 0.0_wp
        ! g%rhs1(1, 5) = g%rhs1(1, 2); g%rhs1(1, 2) = 0.0_wp
        ! g%rhs1(imax, 1) = g%rhs1(imax, 4); g%rhs1(imax, 4) = 0.0_wp
        ! g%nb_diag_1 = [3, 5]

        ! In FDM_Initialize, we do LU decomposition and I need to call this again
        print *, new_line('a'), '1. order, Jacobian 6'
        call FDM_C1N6_Jacobian(imax, g%jac, g%lu1(:, :), g%rhs1(:, :), coef, g%periodic)
        g%nb_diag_1 = [3, 5]

        do ip = 1, size(cases_new)
            ibc = cases_new(ip)
            print *, new_line('a'), 'Bcs case ', ibc

            call FDM_Int1_Initialize(ibc, g%lu1(:, 1:g%nb_diag_1(1)), g%rhs1(:, 1:g%nb_diag_1(2)), lambda, lhs_int, rhs_int)

            nmin = 1
            nmax = imax
            select case (ibc)
            case (BCS_MIN)
                nmin = nmin + 1
                w_n(:, 1) = u(:, 1)
            case (BCS_MAX)
                nmax = nmax - 1
                w_n(:, imax) = u(:, imax)
            end select
            nsize = nmax - nmin + 1

            ! LU decomposition
            select case (g%nb_diag_1(2))
            case (3)
                call TRIDFS(nsize, lhs_int(nmin:nmax, 1), lhs_int(nmin:nmax, 2), lhs_int(nmin:nmax, 3))
            case (5)
               call PENTADFS(nsize, lhs_int(nmin:nmax, 1), lhs_int(nmin:nmax, 2), lhs_int(nmin:nmax, 3), lhs_int(nmin:nmax, 4), lhs_int(nmin:nmax, 5))
            end select

            ! Particular solution
            select case (g%nb_diag_1(1))
            case (3)
                call MatMul_3d(nsize, len, rhs_int(nmin:nmax, 1), rhs_int(nmin:nmax, 3), f(:, nmin:nmax), w_n(:, nmin:nmax))
            case (5)

            end select

            ! BC corrections; to be put inside of new version of matmul_3d
            idr = g%nb_diag_1(2)/2 + 1
            select case (ibc)
            case (BCS_MIN)                    ! BCs at the bottom
                do i = 1, idr - 1
                    w_n(:, 1 + i) = w_n(:, 1 + i) + lhs_int(1 + i, idr - i)*w_n(:, 1)
                end do
            case (BCS_MAX)                    ! BCs at the top
                do i = 1, idr - 1
                    w_n(:, imax - i) = w_n(:, imax - i) + lhs_int(imax - i, idr + i)*w_n(:, imax)
                end do
            end select

            select case (g%nb_diag_1(2))
            case (3)
                call TRIDSS(nsize, len, lhs_int(nmin:nmax, 1), lhs_int(nmin:nmax, 2), lhs_int(nmin:nmax, 3), w_n(:, nmin:nmax))
            case (5)
                call PENTADSS(nsize, len, lhs_int(nmin:nmax, 1), lhs_int(nmin:nmax, 2), lhs_int(nmin:nmax, 3), lhs_int(nmin:nmax, 4), lhs_int(nmin:nmax, 5), w_n(:, nmin:nmax))
            end select

            call check(u, w_n, 'integral.dat')

        end do

! ###################################################################
! Second order equation
! ###################################################################
! to be rewritten in terms of INT_C2NX_
    case (2)

        ! call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g, u, f, tmp)
        ! dummy = lambda*lambda
        ! f = du2_a - dummy*u

        ! ! call INT_C2N6_LHS_E(imax, g%jac, ibc, dummy, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
        ! call PENTADFS(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
        ! call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
        ! call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
        ! ! call INT_C2N6_RHS(imax, i1, g%jac, f, w_n)
        ! call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), w_n(2))
        ! w_n(:) = w_n(:) + u(1)*wrk1d(:, 6) + u(imax)*wrk1d(:, 7)

    end select

    stop

    ! ###################################################################
contains
    subroutine check(u, w_n, name)
        real(wp), intent(in) :: u(:, :), w_n(:, :)
        character(len=*), optional :: name

        real(wp) dummy, error_l2, error_max

        if (present(name)) then
            open (20, file=name)
        end if
        error_l2 = 0.0_wp
        error_max = 0.0_wp
        dummy = 0.0_wp
        do i = 1, size(u, 2)
            do l = 1, size(u, 1)
                if (present(name)) then
                    write (20, 1000) g%nodes(i), u(l, i), w_n(l, i), u(l, i) - w_n(l, i)
                end if
                dummy = dummy + u(l, i)*u(l, i)
                error_l2 = error_l2 + (u(l, i) - w_n(l, i))**2.0_wp
                error_max = max(error_max, abs(u(l, i) - w_n(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(g%jac(1, 1)*dummy)/real(len, wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error_l2)/maxval(abs(u))
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(u))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program VINTEGRAL