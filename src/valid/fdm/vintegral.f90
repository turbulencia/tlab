#include "dns_const.h"

program VINTEGRAL
    use TLab_Constants, only: wp, wi, fmt_r
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLab_Arrays, only: wrk1d, wrk2d, txc
    use FDM, only: fdm_dt, FDM_Initialize
    use FDM, only: FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM4_DIRECT, FDM_COM6_DIRECT
    use FDM_Integral
    use OPR_PARTIAL
    use OPR_ODES

    implicit none

    integer(wi) :: i, l, len

    real(wp), dimension(:, :), pointer :: u => null(), w_n => null(), f => null()
    real(wp), dimension(:, :), pointer :: du1_a => null(), du2_a => null(), du1_n => null(), du2_n => null(), dw1_n => null(), dw2_n => null()
    integer(wi) bcs_aux(2, 2)
    real(wp) :: lambda, wk, x_0
    integer(wi) :: test_type, ibc, ib, im, ndr, ndl

    integer :: bcs_cases(4), fdm_cases(5)
    real(wp), allocatable :: bcs(:, :), x(:), si(:, :)
    character(len=32) :: fdm_names(5)

    type(fdm_dt) :: g
    type(fdm_integral_dt) :: fdmi(2), fdmi_test(2), fdmi_test_lambda(2)

! ###################################################################
! Initialize
    imax = 1
    jmax = 3
    kmax = 256
    len = imax*jmax

    g%size = kmax
    g%scale = 1.0_wp

    isize_field = imax*jmax*kmax
    isize_txc_field = isize_field
    isize_wrk3d = isize_txc_field
    isize_wrk1d = kmax
    isize_wrk2d = max(imax*jmax, max(imax*kmax, jmax*kmax))
    inb_wrk1d = 20
    inb_wrk2d = 3
    inb_txc = 9

    call TLab_Initialize_Memory(__FILE__)
    allocate (x(kmax))

    u(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 1)
    w_n(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 2)
    f(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 3)

    du1_a(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    du2_a(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    du1_n(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    du2_n(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    dw1_n(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 8)
    dw2_n(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 9)

    print *, '1. First order equation.'
    print *, '2. Second order equation; factorized, singular.'
    print *, '3. Second order equation; factorized, regular.'
    print *, '4. Second order equation; direct.'
    print *, '5. Array properties.'
    read (*, *) test_type

    ! ###################################################################
    if (g%periodic) then
        do i = 1, kmax
            x(i) = real(i - 1, wp)/real(kmax, wp)*g%scale
        end do
    else
        ! do i = 1, kmax
        !     x(i) = real(i - 1, wp)/real(kmax - 1, wp)*g%scale
        ! end do
        open (21, file='y.dat')
        do i = 1, kmax
            read (21, *) x(i)
        end do
        close (21)
    end if

    g%mode_fdm1 = FDM_COM6_JACOBIAN     ! default
    g%mode_fdm2 = g%mode_fdm1
    call FDM_Initialize(x, g, fdmi)
    ndr = g%nb_diag_1(2)
    ndl = g%nb_diag_1(1)

    bcs_aux = 0

! ###################################################################
! Define the function f and analytic derivatives
    x_0 = 0.75_wp
    wk = 1.0_wp

    do i = 1, kmax
! single-mode
        ! u(:, i) = 1.0_wp + sin(2.0_wp*pi_wp/g%scale*wk*g%nodes(i)) ! + pi_wp/4.0_wp)
        ! du1_a(:, i) = (2.0_wp*pi_wp/g%scale*wk) &
        !               *cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))! + pi_wp/4.0_wp)
! Gaussian
        u(:, i) = exp(-(g%nodes(i) - x_0*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        du1_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*u(:, i)
        du2_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*du1_a(:, i) &
                      - 1.0_wp/(g%scale/wk)**2*u(:, i)
! exponential
        ! u(:, i) = exp(-g%nodes(i)*wk)
        ! du1_a(:, i) = -wk*u(:, i)
! step
        ! u(:, i) = max(0.0_wp, (g%nodes(i) - g%nodes(kmax/2))*x_0)
        ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, g%nodes(i) - g%nodes(kmax/2)))*0.5_wp*x_0
! tanh
        ! u(:, i) = x_0*log(1.0_wp + exp((g%nodes(i) - g%nodes(kmax/2))/x_0))
        ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(g%nodes(i) - g%nodes(kmax/2))/x_0))
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
    ! We reinitialize integral operator for different FDM schemes, not only the default one
    ! ###################################################################
    select case (test_type)
    case (1)
        write (*, *) 'Eigenvalue ?'
        read (*, *) lambda

        fdm_cases(1:5) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM4_DIRECT, FDM_COM6_DIRECT]
        im = 0
        im = im + 1; fdm_names(im) = 'Jacobian 4'
        im = im + 1; fdm_names(im) = 'Jacobian 6'
        im = im + 1; fdm_names(im) = 'Jacobian 6, penta-diagonal'
        im = im + 1; fdm_names(im) = 'Direct 4'
        im = im + 1; fdm_names(im) = 'Direct 6'
        bcs_cases(1:2) = [BCS_MIN, BCS_MAX]

        do im = 1, 5 !size(fdm_cases)
            print *, new_line('a'), fdm_names(im)

            g%mode_fdm1 = fdm_cases(im)
            g%mode_fdm2 = FDM_COM4_JACOBIAN     ! not used
            call FDM_Initialize(x, g)
            ndr = g%nb_diag_1(2)
            ndl = g%nb_diag_1(1)

            f = du1_a
            ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, f)
            f = f + lambda*u

            do ib = 1, 2
                ibc = bcs_cases(ib)
                print *, new_line('a'), 'Bcs case ', ibc

                fdmi(ib)%bc = ibc
                fdmi(ib)%mode_fdm = g%mode_fdm1
                call FDM_Int1_Initialize(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), lambda, fdmi(ib))

                ! bcs
                select case (ibc)
                case (BCS_MIN)
                    w_n(:, 1) = u(:, 1)
                case (BCS_MAX)
                    w_n(:, kmax) = u(:, kmax)
                end select

                call FDM_Int1_Solve(len, fdmi(ib), fdmi(ib)%rhs, f, w_n, wrk2d, dw1_n(:, 1))

                call check(u, w_n, 'integral.dat')

                ! check the calculation of the derivative at the boundary
                print *, dw1_n(:, 1)
                call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, w_n, dw1_n)
                select case (ibc)
                case (BCS_MIN)
                    print *, dw1_n(:, 1)
                case (BCS_MAX)
                    print *, dw1_n(:, kmax)
                end select

            end do

        end do

! ###################################################################
! Second order equation; singular cases
! ###################################################################
    case (2)
        allocate (bcs(len, 2))
        do i = kmax, 1, -1     ! set the lower value to zero, which is assumed in BCS_NN
            u(:, i) = u(:, i) - u(:, 1)
        end do

        ! f = du2_a
        ! du1_n = du1_a ! I need it for the boundary conditions
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, du1_n)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, du1_n, du2_n)
        ! call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_aux, g, u, du2_n, du1_n)
        f = du2_n

        bcs_cases(1:4) = [BCS_DD, BCS_DN, BCS_ND, BCS_NN]

        do ib = 1, 4
            ibc = bcs_cases(ib)
            print *, new_line('a')

            select case (ibc)
            case (BCS_DD)
                print *, 'Dirichlet/Dirichlet'
                bcs(:, 1) = u(:, 1); bcs(:, 2) = u(:, kmax)
                ! call OPR_ODE2_Factorize_1_SINGULAR_DD_OLD(g%mode_fdm1, g%size, len, g%nodes, g%jac, w_n, f, bcs, dw1_n, wrk1d)
                call OPR_ODE2_Factorize_DD_Sing(len, fdmi, w_n, f, bcs, dw1_n, wrk1d, wrk2d)
            case (BCS_DN)
                print *, 'Dirichlet/Neumann'
                bcs(:, 1) = u(:, 1); bcs(:, 2) = du1_n(:, kmax)
                ! call OPR_ODE2_Factorize_1_SINGULAR_DN_OLD(g%mode_fdm1, g%size, len, g%jac, w_n, f, bcs, dw1_n, wrk1d)
                call OPR_ODE2_Factorize_DN_Sing(len, fdmi, w_n, f, bcs, dw1_n, wrk1d, wrk2d)
            case (BCS_ND)
                print *, 'Neumann/Dirichlet'
                bcs(:, 1) = du1_n(:, 1); bcs(:, 2) = u(:, kmax)
                ! call OPR_ODE2_Factorize_1_SINGULAR_ND_OLD(g%mode_fdm1, g%size, len, g%jac, w_n, f, bcs, dw1_n, wrk1d)
                call OPR_ODE2_Factorize_ND_Sing(len, fdmi, w_n, f, bcs, dw1_n, wrk1d, wrk2d)
            case (BCS_NN)
                print *, 'Neumann/Neumann'
                bcs(:, 1) = du1_n(:, 1); bcs(:, 2) = du1_n(:, kmax)
                ! call OPR_ODE2_Factorize_1_SINGULAR_NN_OLD(g%mode_fdm1, g%size, len, g%jac, w_n, f, bcs, dw1_n, wrk1d)
                call OPR_ODE2_Factorize_NN_Sing(len, fdmi, w_n, f, bcs, dw1_n, wrk1d, wrk2d)
            end select

            call check(u, w_n, 'integral.dat')
            call check(du1_n, dw1_n)

        end do

! ###################################################################
! Second order equation; regular cases
! ###################################################################
    case (3)
        write (*, *) 'Eigenvalue ?'
        read (*, *) lambda

        fdmi(BCS_MIN)%bc = BCS_MIN
        fdmi(BCS_MIN)%mode_fdm = g%mode_fdm1
        call FDM_Int1_Initialize(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), lambda, fdmi(BCS_MIN))

        fdmi(BCS_MAX)%bc = BCS_MAX
        fdmi(BCS_MAX)%mode_fdm = g%mode_fdm1
        call FDM_Int1_Initialize(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), -lambda, fdmi(BCS_MAX))

        allocate (bcs(len, 2))
        ! call random_seed()
        ! call random_number(u)

        ! f = du2_a - lambda*lambda*u
        ! du1_n = du1_a ! I need it for the boundary conditions
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, du1_n)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, du1_n, du2_n)
        ! call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_aux, g, u, du2_n, du1_n)
        f = du2_n - lambda*lambda*u

        bcs_cases(1:2) = [BCS_DD, BCS_NN]!, BCS_DN, BCS_ND]

        do ib = 1, 2
            ibc = bcs_cases(ib)
            print *, new_line('a')

            select case (ibc)
            case (BCS_DD)
                print *, 'Dirichlet/Dirichlet'
                bcs(:, 1) = u(:, 1); bcs(:, 2) = u(:, kmax)
                ! call OPR_ODE2_Factorize_1_REGULAR_DD_OLD(g%mode_fdm1, g%size, len, lambda*lambda, g%jac, w_n, f, bcs, dw1_n, wrk1d)
                call OPR_ODE2_Factorize_DD(len, fdmi, fdmi(BCS_MIN)%rhs, fdmi(BCS_MAX)%rhs, w_n, f, bcs, dw1_n, wrk1d, wrk2d)
            case (BCS_DN) ! not yet developed
                print *, 'Dirichlet/Neumann'
                bcs(:, 1) = u(:, 1); bcs(:, 2) = du1_n(:, kmax)
            case (BCS_ND)
                print *, 'Neumann/Dirichlet'
                bcs(:, 1) = du1_n(:, 1); bcs(:, 2) = u(:, kmax)
            case (BCS_NN)
                print *, 'Neumann/Neumann'
                bcs(:, 1) = du1_n(:, 1); bcs(:, 2) = du1_n(:, kmax)
                ! call OPR_ODE2_Factorize_1_REGULAR_NN_OLD(g%mode_fdm1, g%size, len, lambda*lambda, g%jac, w_n, f, bcs, dw1_n, wrk1d)
                call OPR_ODE2_Factorize_NN(len, fdmi, fdmi(BCS_MIN)%rhs, fdmi(BCS_MAX)%rhs, w_n, f, bcs, dw1_n, wrk1d, wrk2d)
            end select

            call check(u, w_n, 'integral.dat')
            call check(du1_n, dw1_n)

        end do

! ###################################################################
! Second order equation; direct
! ###################################################################
    case (4)
        write (*, *) 'Eigenvalue ?'
        read (*, *) lambda

        allocate (si(g%size, 2))
        allocate (bcs(len, 2))

        g%mode_fdm1 = FDM_COM6_JACOBIAN
        g%mode_fdm2 = FDM_COM6_DIRECT
        call FDM_Initialize(x, g)
        ndr = g%nb_diag_2(2)
        ndl = g%nb_diag_2(1)

        ! call random_seed()
        ! call random_number(u)

        ! f = du2_a - lambda*u
        ! du1_n = du1_a ! I need it for the boundary conditions
        ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, du1_n)
        ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_aux, g, du1_n, du2_n)
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_aux, g, u, du2_n, du1_n)
        f = du2_n - lambda*u

        bcs_cases(1:4) = [BCS_DD, BCS_NN, BCS_DN, BCS_ND]

        do ib = 1, 4
            ibc = bcs_cases(ib)
            print *, new_line('a')

            select case (ibc)
            case (BCS_DD)
                print *, 'Dirichlet/Dirichlet'
                w_n(:, 1) = u(:, 1); w_n(:, kmax) = u(:, kmax)
            case (BCS_DN)
                print *, 'Dirichlet/Neumann'
                w_n(:, 1) = u(:, 1); w_n(:, kmax) = du1_n(:, kmax)
            case (BCS_ND)
                print *, 'Neumann/Dirichlet'
                w_n(:, 1) = du1_n(:, 1); w_n(:, kmax) = u(:, kmax)
            case (BCS_NN)
                print *, 'Neumann/Neumann'
                w_n(:, 1) = du1_n(:, 1); w_n(:, kmax) = du1_n(:, kmax)
            end select

            fdmi(2)%bc = ibc
            call FDM_Int2_Initialize(g%nodes(:), g%lhs2(:, 1:ndl), g%rhs2(:, 1:ndr), lambda, fdmi(2))

            call FDM_Int2_Solve(len, fdmi(2), fdmi(2)%rhs, f, w_n, wrk2d)
            call check(u, w_n, 'integral.dat')

        end do

! ###################################################################
! Test properties of integral matrices
! ###################################################################
    case (5)
        write (*, *) 'Eigenvalue ?'
        read (*, *) lambda

        fdm_cases(1:5) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM4_DIRECT, FDM_COM6_DIRECT]
        fdm_names(1:2) = ['1. order, Jacobian 4', '1. order, Jacobian 6']
        fdm_names(3) = '1. order, Jacobian 6 Penta'
        fdm_names(4:5) = ['1. order, Direct 4', '1. order, Direct 6']
        bcs_cases(1:2) = [BCS_MIN, BCS_MAX]

        do im = 1, 5 !size(fdm_cases)
            print *, new_line('a'), fdm_names(im)

            g%mode_fdm1 = fdm_cases(im)
            g%mode_fdm2 = FDM_COM4_JACOBIAN ! not used
            call FDM_Initialize(x, g)
            ndr = g%nb_diag_1(2)
            ndl = g%nb_diag_1(1)

            do ib = 1, 2
                ibc = bcs_cases(ib)
                print *, new_line('a'), 'Bcs case ', ibc

                fdmi(ib)%bc = ibc
                fdmi(ib)%mode_fdm = g%mode_fdm1
                call FDM_Int1_CreateSystem(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), 0.0_wp, fdmi(ib))

                fdmi_test(ib)%bc = fdmi(ib)%bc
                fdmi_test(ib)%mode_fdm = fdmi(ib)%mode_fdm
                call FDM_Int1_CreateSystem(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), 1.0_wp, fdmi_test(ib))

                call check(fdmi(ib)%rhs, fdmi_test(ib)%rhs, 'integral.dat')
                ! print*,fdmi(ib)%rhs_b
                ! print*,fdmi_test(ib)%rhs_b
                ! print*,fdmi(ib)%rhs_t
                ! print*,fdmi_test(ib)%rhs_t

                ! checking linearity in lhs
                fdmi_test_lambda(ib)%bc = fdmi(ib)%bc
                fdmi_test_lambda(ib)%mode_fdm = fdmi(ib)%mode_fdm
                call FDM_Int1_CreateSystem(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), lambda, fdmi_test_lambda(ib))

                call check(fdmi(ib)%lhs + lambda*(fdmi_test(ib)%lhs - fdmi(ib)%lhs), &
                           fdmi_test_lambda(ib)%lhs, 'integral.dat')

            end do

        end do

    end select

    stop

    ! ###################################################################
contains
    subroutine check(u, w_n, name)
        real(wp), intent(in) :: u(:, :), w_n(:, :)
        character(len=*), optional :: name

        real(wp) dummy, error_l2, error_max
        character(len=*), parameter :: fmt_loc = '(5(1x, '//fmt_r//'))'

        if (present(name)) then
            open (20, file=name)
        end if
        error_l2 = 0.0_wp
        error_max = 0.0_wp
        dummy = 0.0_wp
        do i = 1, size(u, 2)
            do l = 1, size(u, 1)
                if (present(name)) then
                    write (20, fmt_loc) g%nodes(i), u(l, i), w_n(l, i), u(l, i) - w_n(l, i)
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
    end subroutine check

end program VINTEGRAL
