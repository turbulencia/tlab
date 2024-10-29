#include "dns_const.h"

program VPARTIAL
    use TLab_Constants
    use TLab_Types, only: grid_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLAB_VARS, only: visc, schmidt
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLab_Arrays, only: wrk1d, wrk2d, txc, x!, wrk3d
    use FDM_ComX_Direct
    use FDM_PROCS
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use OPR_PARTIAL

    implicit none

    type(grid_dt) :: g

    integer(wi) :: i, l, len

    real(wp), dimension(:, :), pointer :: u
    real(wp), dimension(:, :), pointer :: du1_a, du1_b, du1_c, du1_n
    real(wp), dimension(:, :), pointer :: du2_a, du2_n1, du2_n2, du2_n3
    integer(wi) bcs_aux(2, 2)
    real(wp) :: wk, coef(5), dummy
    integer(wi) :: test_type, ibc, ip, ic, ndr, idr, ndl, idl, im
    integer(wi) :: nmin, nmax, nsize

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
    du1_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 2)
    du1_b(1:len, 1:imax) => txc(1:imax*jmax*kmax, 3)
    du1_c(1:len, 1:imax) => txc(1:imax*jmax*kmax, 4)
    du1_n(1:len, 1:imax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:len, 1:imax) => txc(1:imax*jmax*kmax, 6)
    du2_n1(1:len, 1:imax) => txc(1:imax*jmax*kmax, 7)
    du2_n2(1:len, 1:imax) => txc(1:imax*jmax*kmax, 8)
    du2_n3(1:len, 1:imax) => txc(1:imax*jmax*kmax, 9)

    call TLab_Allocate_Real(__FILE__, x, [g%size, g%inb_grid], g%name)

    ! Valid settings
    test_type = 1

    g%periodic = .false.
    ! g%periodic = .true.
    wk = 1.0_wp ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) wk
    g%mode_fdm1 = FDM_COM6_JACOBIAN
    ! g%mode_fdm1 = FDM_COM6_JACOBIAN_PENTA
    ! g%mode_fdm1 = FDM_COM6_DIRECT
    g%mode_fdm2 = g%mode_fdm1

!  ###################################################################

    if (g%periodic) then
        do i = 1, imax
            x(i, 1) = real(i - 1, wp)/real(imax, wp)*g%scale
        end do
    else
        do i = 1, imax
            x(i, 1) = real(i - 1, wp)/real(imax - 1, wp)*g%scale
        end do
        ! open (21, file='y.dat')
        ! do i = 1, imax
        !     read (21, *) x(i, 1)
        ! end do
        ! close (21)
        ! g%scale = x(imax, 1) - x(1, 1)
    end if

    call FDM_INITIALIZE(x, g, wrk1d)

! Bcs
    bcs_aux = 0

! ###################################################################
! Define the function and analytic derivatives
    do i = 1, imax
        do l = 1, len
! single-mode
            ! u(l, i) = 1.0_wp + &
            !           sin(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))!+pi_wp/C_4_R)
            ! du1_a(l, i) = (2.0_wp*pi_wp/g%scale*wk) &
            !               *cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))!+pi_wp/C_4_R)
            u(l, i) = 1.0_wp + &
                      cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))!+pi_wp/C_4_R)
            du1_a(l, i) = -(2.0_wp*pi_wp/g%scale*wk) &
                          *sin(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))!+pi_wp/C_4_R)

            du2_a(l, i) = -(2.0_wp*pi_wp/g%scale*wk)**2 &
                          *cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))!+pi_wp/C_4_R)

! Gaussian
            ! dummy = 1.0_wp / ( 2.0_wp*(g%scale/M_REAL(wk*l))**2 )
            ! u(l,i)     = EXP(-dummy*(g%nodes(i)-x_0*g%scale)**2)
            ! du1_a(l,i) =-2.0_wp *dummy *(g%nodes(i)-x_0*g%scale) *u(l,i)
            ! du2_a(l,i) =-2.0_wp *dummy *(g%nodes(i)-x_0*g%scale) *du1_a(l,i) - 2.0_wp *dummy *u(l,i)
! Exponential
            ! u(l,i)     = EXP(g%nodes(i)/(g%scale/wk))
            ! du1_a(l,i) = wk/g%scale*u(l,i)
            ! du2_a(l,i) = wk/g%scale*du1_a(l,i)
! delta-function
            ! u(i)     = MAX(0.0_wp,2.0_wp-M_REAL(i))
            ! du1_a(i) = 0.0_wp
            ! du2_a(i) = 0.0_wp
! hyperboic tangent
            ! u(l,i)     = wk*LOG(1.0_wp+EXP(g%nodes(i)/wk))
            ! du1_a(l,i) = C_05_R*(1.0_wp+TANH(C_05_R*g%nodes(i)/wk))
            ! du2_a(l,i) = C_025_R/wk/(COSH(C_05_R*g%nodes(i)/wk))**2
! Polynomial
            ! dummy = 4.0_wp
            ! u(l, i) = ((g%scale - g%nodes(i))/wk)**dummy
            ! du1_a(l, i) = -dummy/wk*((g%scale - g%nodes(i))/wk)**(dummy - 1.0_wp)
            ! du2_a(l, i) = dummy*(dummy - 1.0_wp)*((g%scale - g%nodes(i))/wk)**(dummy - 2.0_wp)
        end do
    end do

! ###################################################################
! First-order derivative
! ###################################################################
    select case (test_type)
    case (1)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, du1_b)

        ! Jacobian based
        print *, '1. order, Jacobian 4'
        call FDM_C1N4_Jacobian(g%size, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)
        call TRIDFS(g%size, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3))
        call MatMul_3d_antisym(g%size, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), u, du1_n, g%periodic)
        call TRIDSS(g%size, len, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), du1_n)
        call check(u, du1_a, du1_n, 'partial.dat')

        print *, '1. order, Jacobian 6'
        call FDM_C1N6_Jacobian(g%size, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)
        call TRIDFS(g%size, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3))
        call MatMul_5d_antisym(g%size, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), u, du1_n, g%periodic)
        call TRIDSS(g%size, len, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), du1_n)
        ! call TRIDPFS(g%size, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), g%lu1(:, 4), g%lu1(:, 5))
        ! call TRIDPSS(g%size, len, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), g%lu1(:, 4), g%lu1(:, 5), du1_n, wrk2d)
        call check(u, du1_a, du1_n, 'partial.dat')

        print *, '1. order, Jacobian 6 Penta'
        call FDM_C1N6_Jacobian_Penta(g%size, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)
        call PENTADFS2(g%size, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), g%lu1(:, 4), g%lu1(:, 5))
     call MatMul_7d_antisym(g%size, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), u, du1_n, g%periodic)
        call PENTADSS2(g%size, len, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), g%lu1(:, 4), g%lu1(:, 5), du1_n)
        ! call PENTADPFS(g%size, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), g%lu1(:, 4), g%lu1(:, 5), g%lu1(:, 6), g%lu1(:, 7))
        ! call PENTADPSS(g%size, len, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), g%lu1(:, 4), g%lu1(:, 5), g%lu1(:, 6), g%lu1(:, 7), du1_n)
        call check(u, du1_a, du1_n, 'partial.dat')

        ! Direct metrics
        call FDM_C1N4_Direct(g%size, x, g%lu1, g%rhs1, g%nb_diag_1)

        call FDM_C1N6_Direct(g%size, x, g%lu1, g%rhs1, g%nb_diag_1)

        ! -------------------------------------------------------------------
        !   Testing the reduction routines
        fdm_cases(1:3) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA]
        fdm_names(1:2) = ['1. order, Jacobian 4', '1. order, Jacobian 6']; fdm_names(3) = '1. order, Jacobian 6 Penta'
        bcs_cases(1:3) = [BCS_MIN, BCS_MAX, BCS_BOTH]
        do ip = 1, 3
            ibc = bcs_cases(ip)
            print *, new_line('a'), 'Bcs case ', ibc

            nmin = 1
            nmax = g%size
            if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                nmin = nmin + 1
            end if
            if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                nmax = nmax - 1
            end if
            nsize = nmax - nmin + 1

            do im = 1, 3
                g%mode_fdm1 = fdm_cases(im)
                print *, fdm_names(im)

                select case (g%mode_fdm1)
                case (FDM_COM4_JACOBIAN)
                    call FDM_C1N4_Jacobian(imax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN)
                    call FDM_C1N6_Jacobian(imax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN_PENTA)
                    call FDM_C1N6_Jacobian_Penta(imax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                end select
                ndl = g%nb_diag_1(1)
                idl = g%nb_diag_1(1)/2 + 1
                ndr = g%nb_diag_1(2)
                idr = g%nb_diag_1(2)/2 + 1

                ! g%rhsr_b = 0.0_wp
                ! g%rhsr_t = 0.0_wp
                call FDM_Bcs_Reduce(ibc, g%lu1(:, 1:ndl), g%rhs1(:, 1:ndr), g%rhsr_b, g%rhsr_t)

                select case (g%nb_diag_1(1))
                case (3)
                    call TRIDFS(nsize, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3))
                case (5)
                    call PENTADFS2(nsize, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), g%lu1(nmin:nmax, 4), g%lu1(nmin:nmax, 5))
                end select

                du1_n(:, 1) = u(:, 1)           ! boundary condition
                du1_n(:, imax) = u(:, imax)
                select case (g%nb_diag_1(2))
                case (3)
                    call MatMul_3d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=g%rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                case (5)
                    call MatMul_5d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=g%rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                case (7)
                    call MatMul_7d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=g%rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                end select

                select case (g%nb_diag_1(1))
                case (3)
                    call TRIDSS(nsize, len, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), du1_n(:, nmin:nmax))
                case (5)
                    call PENTADSS2(nsize, len, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), g%lu1(nmin:nmax, 4), g%lu1(nmin:nmax, 5), du1_n(:, nmin:nmax))
                end select

                if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                    du1_n(:, 1) = wrk2d(:, 1)
                    do ic = 1, idl - 1
                        du1_n(:, 1) = du1_n(:, 1) + g%lu1(1, idl + ic)*du1_n(:, 1 + ic)
                    end do
                    du1_n(:, 1) = du1_n(:, 1) + g%lu1(1, 1)*du1_n(:, 1 + ic)
                end if
                if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                    du1_n(:, imax) = wrk2d(:, 2)
                    do ic = 1, idl - 1
                        du1_n(:, imax) = du1_n(:, imax) + g%lu1(imax, idl - ic)*du1_n(:, imax - ic)
                    end do
                    du1_n(:, imax) = du1_n(:, imax) + g%lu1(imax, ndl)*du1_n(:, imax - ic)
                end if
                call check(u, du1_a, du1_n, 'partial.dat')

            end do

        end do

! ###################################################################
! Boundary conditions
! ###################################################################
    case (2)
        fdm_cases(1:3) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA]
        fdm_names(1:2) = ['1. order, Jacobian 4', '1. order, Jacobian 6']; fdm_names(3) = '1. order, Jacobian 6 Penta'
        bcs_cases(1:4) = [BCS_DD, BCS_ND, BCS_DN, BCS_NN]
        do ip = 1, 4
            ibc = bcs_cases(ip)
            print *, new_line('a'), 'Bcs case ', ibc

            nmin = 1
            nmax = g%size
            if (any([BCS_ND, BCS_NN] == ibc)) then
                nmin = nmin + 1
            end if
            if (any([BCS_DN, BCS_NN] == ibc)) then
                nmax = nmax - 1
            end if
            nsize = nmax - nmin + 1

            du1_n(:, 1) = du1_a(:, 1)
            du1_n(:, imax) = du1_a(:, imax)

            do im = 1, 3
                g%mode_fdm1 = fdm_cases(im)
                print *, fdm_names(im)

                select case (g%mode_fdm1)
                case (FDM_COM4_JACOBIAN)
                    call FDM_C1N4_Jacobian(imax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN)
                    call FDM_C1N6_Jacobian(imax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN_PENTA)
                    call FDM_C1N6_Jacobian_Penta(imax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                end select
                idl = g%nb_diag_1(1)/2 + 1
                idr = g%nb_diag_1(2)/2 + 1

                call FDM_Bcs_Neumann(ibc, g%lu1(:, 1:g%nb_diag_1(1)), g%rhs1(:, 1:g%nb_diag_1(2)), g%rhs1_b, g%rhs1_t)

                select case (g%nb_diag_1(1))
                case (3)
                    call TRIDFS(nsize, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3))
                case (5)
                    call PENTADFS2(nsize, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), g%lu1(nmin:nmax, 4), g%lu1(nmin:nmax, 5))
                end select

                select case (g%nb_diag_1(2))
                case (3)
                    call MatMul_3d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhs1_b, bcs_b=wrk2d(:, 1), rhs_t=g%rhs1_t, bcs_t=wrk2d(:, 2))
                case (5)
                    call MatMul_5d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhs1_b, bcs_b=wrk2d(:, 1), rhs_t=g%rhs1_t, bcs_t=wrk2d(:, 2))
                case (7)
                    call MatMul_7d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhs1_b, bcs_b=wrk2d(:, 1), rhs_t=g%rhs1_t, bcs_t=wrk2d(:, 2))
                end select

                select case (g%nb_diag_1(1))
                case (3)
                    call TRIDSS(nsize, len, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), du1_n(:, nmin:nmax))
                case (5)
                    call PENTADSS2(nsize, len, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), g%lu1(nmin:nmax, 4), g%lu1(nmin:nmax, 5), du1_n(:, nmin:nmax))
                end select

                call check(u(:, nmin:nmax), du1_a(:, nmin:nmax), du1_n(:, nmin:nmax), 'partial.dat')
                if (any([BCS_ND, BCS_NN] == ibc)) then
                    do ic = 1, idl - 1
                        wrk2d(1:len, 1) = wrk2d(1:len, 1) + g%lu1(1, idl + ic)*du1_n(:, 1 + ic)
                    end do
                    print *, u(:, 1), wrk2d(1:len, 1)
                end if
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    do ic = 1, idl - 1
                        wrk2d(1:len, 2) = wrk2d(1:len, 2) + g%lu1(imax, idl - ic)*du1_n(:, imax - ic)
                    end do
                    print *, u(:, imax), wrk2d(1:len, 2)
                end if

            end do

        end do

! ###################################################################
! Second-order derivative
! ###################################################################
    case (3)
        ! call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_aux, g, u, du2_n2, du1_n)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, du1_n, du2_n1)

        ! Jacobian based
        fdm_cases(1:3) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER]
        fdm_names(1:2) = ['2. order, Jacobian 4', '2. order, Jacobian 6']
        fdm_names(3) = '2. order, Jacobian 6 hyper'     ! funny, cannot be in previous line because different length of character string

        do im = 1, 3
            g%mode_fdm2 = fdm_cases(im)
            print *, fdm_names(im)

            select case (g%mode_fdm2)
            case (FDM_COM4_JACOBIAN)
                call FDM_C2N4_Jacobian(imax, g%jac, g%lu2, g%rhs2, g%nb_diag_2, coef, g%periodic)

            case (FDM_COM6_JACOBIAN)
                call FDM_C2N6_Jacobian(imax, g%jac, g%lu2, g%rhs2, g%nb_diag_2, coef, g%periodic)

            case (FDM_COM6_JACOBIAN_HYPER)
                call FDM_C2N6_Hyper_Jacobian(imax, g%jac, g%lu2, g%rhs2, g%nb_diag_2, coef, g%periodic)

            end select

            call TRIDFS(g%size, g%lu2(:, 1), g%lu2(:, 2), g%lu2(:, 3))

            select case (g%nb_diag_2(2))
            case (3)
                call MatMul_3d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), u, du1_n, g%periodic, &
                                       ibc, rhs_b=g%rhs1_b, bcs_b=wrk2d(:, 1), rhs_t=g%rhs1_t, bcs_t=wrk2d(:, 2))
            case (5)
                call MatMul_5d_sym(imax, len, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), u, du2_n2, g%periodic)

            case (7)
          call MatMul_7d_sym(imax, len, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), g%rhs2(:, 6), g%rhs2(:, 7), u, du2_n2, g%periodic)

            end select

            ip = g%nb_diag_2(2)      ! add Jacobian correction A_2 dx2 du
            call MatMul_3d_add(imax, len, g%rhs2(:, ip + 1), g%rhs2(:, ip + 2), g%rhs2(:, ip + 3), du1_a, du2_n2)

            call TRIDSS(g%size, len, g%lu2(:, 1), g%lu2(:, 2), g%lu2(:, 3), du2_n2)
            call check(u, du2_a, du2_n2, 'partial.dat')

        end do

        ! Direct metrics
        print *, '2. order, Direct 4'
        call FDM_C2N4_Direct(imax, x, wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_2)
        call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call MatMul_5d(imax, len, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), u, du2_n2)
        call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du2_n2)
        call check(u, du2_a, du2_n2, 'partial.dat')

        print *, '2. order, Direct 6'
        call FDM_C2N6_Direct(imax, x, wrk1d(:, 1), wrk1d(:, 4), g%nb_diag_2)
        call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call MatMul_5d(imax, len, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), u, du2_n2)
        call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du2_n2)
        call check(u, du2_a, du2_n2, 'partial.dat')

    end select

    stop

    ! ###################################################################
contains
    subroutine check(u, du_a, du_n, name)
        real(wp), intent(in) :: u(:, :), du_a(:, :), du_n(:, :)
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
                    write (20, 1000) g%nodes(i), u(l, i), du_a(l, i), du_n(l, i), du_a(l, i) - du_n(l, i)
                end if
                dummy = dummy + du_a(l, i)*du_a(l, i)
                error_l2 = error_l2 + (du_a(l, i) - du_n(l, i))**2.0_wp
                error_max = max(error_max, abs(du_a(l, i) - du_n(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(g%jac(1, 1)*dummy)/real(len, wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error_l2)/maxval(abs(du1_a))
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(du1_a))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program VPARTIAL
