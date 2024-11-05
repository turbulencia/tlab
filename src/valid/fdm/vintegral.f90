#include "dns_const.h"

program VINTEGRAL
    use TLab_Constants
    use FDM, only: grid_dt, x, FDM_INITIALIZE
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLAB_VARS, only: visc, schmidt
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLab_Arrays, only: wrk1d, txc
    ! use FDM, only: x
    use FDM_ComX_Direct
    use FDM_Integrate
    use FDM_MatMul
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use OPR_PARTIAL
    use OPR_ODES

    implicit none

    type(grid_dt) :: g

    integer(wi) :: i, l, len

    real(wp), dimension(:, :), pointer :: u, w_n, f
    real(wp), dimension(:, :), pointer :: du1_a, dw1_n, du2_a
    integer(wi) bcs_aux(2, 2)
    real(wp) :: lambda, coef(5), wk, x_0!, dummy
    integer(wi) :: test_type, ibc, ib, ip, im, idr, ndr, ndl!, ic

    integer, parameter :: i1 = 1
    integer :: bcs_cases(4), fdm_cases(3)
    character(len=32) :: fdm_names(3)

! ###################################################################
! Initialize
    imax = 129
    jmax = 1
    kmax = 1
    len = jmax*kmax

    visc = 1.0_wp   ! Needed in FDM_INITIALIZE
    schmidt = 1.0_wp

    g%inb_grid = 99
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

    call TLab_Initialize_Memory(__FILE__)

    u(1:len, 1:imax) => txc(1:imax*jmax*kmax, 1)
    w_n(1:len, 1:imax) => txc(1:imax*jmax*kmax, 2)
    f(1:len, 1:imax) => txc(1:imax*jmax*kmax, 3)

    du1_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 4)
    dw1_n(1:len, 1:imax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 6)

    ! call TLab_Allocate_Real(__FILE__, x, [g%size, g%inb_grid], g%name)

    g%periodic = .false.

    wk = 1.0_wp ! WRITE(*,*) 'Wavenumber ?'; READ(*,*) wk
    write (*, *) 'Eigenvalue ?'
    read (*, *) lambda

    test_type = 1

    ! ###################################################################
    if (g%periodic) then
        do i = 1, imax
            ! x(i, 1) = real(i - 1, wp)/real(imax, wp)*g%scale
            wrk1d(i,1) = real(i - 1, wp)/real(imax, wp)*g%scale
        end do
    else
        ! do i = 1, imax
        !     x(i, 1) = real(i - 1, wp)/real(imax - 1, wp)*g%scale
        ! end do
        open (21, file='y.dat')
        do i = 1, imax
            read (21, *) wrk1d(i,1) !x(i, 1)
        end do
        close (21)
        ! g%scale = x(imax, 1) - x(1, 1)
    end if

    ! to calculate the Jacobians
    g%mode_fdm1 = FDM_COM6_JACOBIAN ! FDM_COM6_JACOBIAN_PENTA
    g%mode_fdm2 = g%mode_fdm1
    call FDM_INITIALIZE(x, g, wrk1d, wrk1d(:,4))

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
        fdm_cases(1:3) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA]
        fdm_names(1:2) = ['1. order, Jacobian 4', '1. order, Jacobian 6']; fdm_names(3) = '1. order, Jacobian 6 Penta'
        bcs_cases(1:2) = [BCS_MIN, BCS_MAX]

        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, f)
        f = du1_a + lambda*u

        do ib = 1, 2
            ibc = bcs_cases(ib)
            print *, new_line('a'), 'Bcs case ', ibc

            do im = 1, 3
                g%mode_fdm1 = fdm_cases(im)
                print *, fdm_names(im)

                select case (g%mode_fdm1)
                case (FDM_COM4_JACOBIAN)
                    call FDM_C1N4_Jacobian(imax, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN)
                    call FDM_C1N6_Jacobian(imax, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN_PENTA)
                    call FDM_C1N6_Jacobian_Penta(imax, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                end select
                ! idl = g%nb_diag_1(1)/2 + 1
                idr = g%nb_diag_1(2)/2 + 1
                ndr = g%nb_diag_1(2)
                ndl = g%nb_diag_1(1)

                ip = (ib - 1)*g%nb_diag_1(2)

                call FDM_Int1_Initialize(ibc, g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), lambda, &
                                         g%lhsi(:, ip + 1:ip + ndr), g%rhsi(:, (ib - 1)*ndl + 1:(ib - 1)*ndl + ndl), &
                                         g%rhsi_b((ib - 1)*5 + 1:, :), g%rhsi_t((ib - 1)*5:, :))
                ! LU decomposition
                select case (g%nb_diag_1(2))
                case (3)
                    call TRIDFS(g%size - 2, g%lhsi(2:, ip + 1), g%lhsi(2:, ip + 2), g%lhsi(2:, ip + 3))
                case (5)
                    call PENTADFS(g%size - 2, g%lhsi(2:, ip + 1), g%lhsi(2:, ip + 2), g%lhsi(2:, ip + 3), &
                                  g%lhsi(2:, ip + 4), g%lhsi(2:, ip + 5))
                case (7)
                    call HEPTADFS(g%size - 2, g%lhsi(2:, ip + 1), g%lhsi(2:, ip + 2), g%lhsi(2:, ip + 3), &
                                  g%lhsi(2:, ip + 4), g%lhsi(2:, ip + 5), g%lhsi(2:, ip + 6), g%lhsi(2:, ip + 7))
                end select

                ! bcs
                select case (ibc)
                case (BCS_MIN)
                    w_n(:, 1) = u(:, 1)
                case (BCS_MAX)
                    w_n(:, imax) = u(:, imax)
                end select

                call OPR_Integral1(len, g, f, w_n, ibc)

                call check(u, w_n, 'integral.dat')

            end do

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
