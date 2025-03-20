program VPARTIAL
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLab_Arrays, only: wrk2d, txc
    use FDM, only: fdm_dt, FDM_Initialize, FDM_Der1_Solve, FDM_Der2_Solve
    use FDM, only: FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM6_JACOBIAN_HYPER, FDM_COM4_DIRECT, FDM_COM6_DIRECT
    use FDM_ComX_Direct
    use FDM_PROCS
    use FDM_MatMul
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian

    implicit none

    integer(wi) :: i, l, len

    real(wp), dimension(:, :), pointer :: u
    real(wp), dimension(:, :), pointer :: du1_a, du1_b, du1_c, du1_n
    real(wp), dimension(:, :), pointer :: du2_a, du2_n1, du2_n2, du2_n3
    integer(wi) bcs_aux(2, 2)
    real(wp) :: wk, x_0, coef(5)!, dummy
    integer(wi) :: test_type, ibc, ip, ic, ndr, idr, ndl, idl, im
    integer(wi) :: nmin, nmax, nsize
    real(wp) rhsr_b(5, 0:7), rhsr_t(0:4, 8)

    integer, parameter :: i1 = 1
    integer :: bcs_cases(4), fdm_cases(5)
    real(wp), allocatable :: x(:)
    character(len=32) :: fdm_names(5)

    type(fdm_dt) g

! ###################################################################
! Initialize
    imax = 2
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

    du1_a(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 2)
    du1_b(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    du1_c(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    du1_n(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    du2_n1(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    du2_n2(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 8)
    du2_n3(1:len, 1:kmax) => txc(1:imax*jmax*kmax, 9)

    print *, '1. First-order derivative.'
    print *, '2. Second-order derivative.'
    print *, '3. Reduction routines.'
    print *, '4. Boundary Conditions.'
    read (*, *) test_type

!  ###################################################################
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
    call FDM_Initialize(x, g)
    ndr = g%nb_diag_1(2)
    ndl = g%nb_diag_1(1)

    bcs_aux = 0

! ###################################################################
! Define the function and analytic derivatives
    x_0 = 0.75_wp
    wk = 1.0_wp

    do i = 1, kmax
        ! ! single-mode
        ! u(:, i) = 1.0_wp + sin(2.0_wp*pi_wp/g%scale*wk*g%nodes(i)) ! + pi_wp/4.0_wp)
        ! du1_a(:, i) = (2.0_wp*pi_wp/g%scale*wk) &
        !               *cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))! + pi_wp/4.0_wp)
        ! Gaussian
        u(:, i) = exp(-(g%nodes(i) - x_0*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        du1_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*u(:, i)
        du2_a(:, i) = -(g%nodes(i) - x_0*g%scale)/(g%scale/wk)**2*du1_a(:, i) &
                      - 1.0_wp/(g%scale/wk)**2*u(:, i)
        ! ! exponential
        ! u(:, i) = exp(-g%nodes(i)*wk)
        ! du1_a(:, i) = -wk*u(:, i)
        ! ! step
        ! u(:, i) = max(0.0_wp, (g%nodes(i) - g%nodes(kmax/2))*x_0)
        ! du1_a(:, i) = (1.0_wp + sign(1.0_wp, g%nodes(i) - g%nodes(kmax/2)))*0.5_wp*x_0
        ! ! tanh
        ! u(:, i) = x_0*log(1.0_wp + exp((g%nodes(i) - g%nodes(kmax/2))/x_0))
        ! du1_a(:, i) = 0.5_wp*(1.0_wp + tanh(0.5_wp*(g%nodes(i) - g%nodes(kmax/2))/x_0))
        ! ! Polynomial
        ! dummy = 4.0_wp
        ! u(:, i) = ((g%scale - g%nodes(i))/wk)**dummy
        ! du1_a(:, i) = -dummy*((g%scale - g%nodes(i))/wk)**(dummy - 1.0_wp)
        ! ! zero
        ! u(:, i) = 0.0_wp
        ! du1_a(:, i) = 0.0_wp
        ! ! delta-function
        ! u(:, i) = max(0.0_wp, 2.0_wp - real(i, wp))
        ! du1_a(:, i) = 0.0_wp
        ! du2_a(:, i) = 0.0_wp
    end do

    ! ###################################################################
    ! First-order derivative
    ! ###################################################################
    select case (test_type)
    case (1)
        fdm_cases(1:5) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM4_DIRECT, FDM_COM6_DIRECT]
        im = 0
        im = im + 1; fdm_names(im) = 'Jacobian 4'
        im = im + 1; fdm_names(im) = 'Jacobian 6'
        im = im + 1; fdm_names(im) = 'Jacobian 6, penta-diagonal'
        im = im + 1; fdm_names(im) = 'Direct 4'     ! undeveloped
        im = im + 1; fdm_names(im) = 'Direct 6'     ! undeveloped

        do im = 1, 4 !size(fdm_cases)
            print *, new_line('a'), fdm_names(im)

            g%mode_fdm1 = fdm_cases(im)
            call FDM_Initialize(x, g)

            call FDM_Der1_Solve(len, bcs_aux(:, 1), g, g%lu1, u, du1_n, wrk2d)

            call check(u, du1_a, du1_n, 'partial.dat')

        end do

        ! ###################################################################
        ! Second-order derivative
        ! ###################################################################
    case (2)
        fdm_cases(1:5) = [FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER, FDM_COM4_DIRECT, FDM_COM6_DIRECT]
        im = 0
        im = im + 1; fdm_names(im) = 'Jacobian 4'
        im = im + 1; fdm_names(im) = 'Jacobian 6'
        im = im + 1; fdm_names(im) = 'Jacobian 6, hyper-diffusive'
        im = im + 1; fdm_names(im) = 'Direct 4'
        im = im + 1; fdm_names(im) = 'Direct 6'

        do im = 1, 5 !size(fdm_cases)
            print *, new_line('a'), fdm_names(im)

            g%mode_fdm1 = fdm_cases(im)
            g%mode_fdm2 = fdm_cases(im)
            if (g%mode_fdm2 == FDM_COM4_DIRECT) g%mode_fdm1 = FDM_COM4_JACOBIAN     ! 1. order are undeveloped
            if (g%mode_fdm2 == FDM_COM6_DIRECT) g%mode_fdm1 = FDM_COM6_JACOBIAN
            call FDM_Initialize(x, g)

            call FDM_Der1_Solve(len, bcs_aux(:, 1), g, g%lu1, u, du1_n, wrk2d)  ! I need du1_n in Jacobian formulation
            call FDM_Der2_Solve(len, g, g%lu2, u, du2_n1, du1_n, wrk2d)

            call check(u, du2_a, du2_n1, 'partial.dat')

        end do

        ! ###################################################################
        !   Testing the reduction routines
        ! ###################################################################
    case (3)
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
                    call FDM_C1N4_Jacobian(kmax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN)
                    call FDM_C1N6_Jacobian(kmax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN_PENTA)
                    call FDM_C1N6_Jacobian_Penta(kmax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                end select
                ndl = g%nb_diag_1(1)
                idl = g%nb_diag_1(1)/2 + 1
                ndr = g%nb_diag_1(2)
                idr = g%nb_diag_1(2)/2 + 1

                ! g%rhsr_b = 0.0_wp
                ! g%rhsr_t = 0.0_wp
                call FDM_Bcs_Reduce(ibc, g%lu1(:, 1:ndl), g%rhs1(:, 1:ndr), rhsr_b, rhsr_t)

                select case (g%nb_diag_1(1))
                case (3)
                    call TRIDFS(nsize, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3))
                case (5)
                    call PENTADFS2(nsize, g%lu1(nmin:nmax, 1), g%lu1(nmin:nmax, 2), g%lu1(nmin:nmax, 3), g%lu1(nmin:nmax, 4), g%lu1(nmin:nmax, 5))
                end select

                du1_n(:, 1) = u(:, 1)           ! boundary condition
                du1_n(:, kmax) = u(:, kmax)
                select case (g%nb_diag_1(2))
                case (3)
                    call MatMul_3d_antisym(kmax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), &
                                           u, du1_n, g%periodic, &
                                           ibc, rhs_b=rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                case (5)
                    call MatMul_5d_antisym(kmax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), &
                                           u, du1_n, g%periodic, &
                                           ibc, rhs_b=rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
                case (7)
                    call MatMul_7d_antisym(kmax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), &
                                           u, du1_n, g%periodic, &
                                           ibc, rhs_b=rhsr_b(:, 1:), bcs_b=wrk2d(:, 1), rhs_t=rhsr_t(1:, :), bcs_t=wrk2d(:, 2))
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
                    du1_n(:, kmax) = wrk2d(:, 2)
                    do ic = 1, idl - 1
                        du1_n(:, kmax) = du1_n(:, kmax) + g%lu1(kmax, idl - ic)*du1_n(:, kmax - ic)
                    end do
                    du1_n(:, kmax) = du1_n(:, kmax) + g%lu1(kmax, ndl)*du1_n(:, kmax - ic)
                end if
                call check(u, du1_a, du1_n, 'partial.dat')

            end do

        end do

        ! ###################################################################
        ! Boundary conditions
        ! ###################################################################
    case (4)
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
            du1_n(:, kmax) = du1_a(:, kmax)

            do im = 1, 3
                g%mode_fdm1 = fdm_cases(im)
                print *, fdm_names(im)

                select case (g%mode_fdm1)
                case (FDM_COM4_JACOBIAN)
                    call FDM_C1N4_Jacobian(kmax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN)
                    call FDM_C1N6_Jacobian(kmax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

                case (FDM_COM6_JACOBIAN_PENTA)
                    call FDM_C1N6_Jacobian_Penta(kmax, g%jac, g%lu1, g%rhs1, g%nb_diag_1, coef, g%periodic)

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
                    call MatMul_3d_antisym(kmax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), &
                                           u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhs1_b, bcs_b=wrk2d(:, 1), rhs_t=g%rhs1_t, bcs_t=wrk2d(:, 2))
                case (5)
                    call MatMul_5d_antisym(kmax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), &
                                           u, du1_n, g%periodic, &
                                           ibc, rhs_b=g%rhs1_b, bcs_b=wrk2d(:, 1), rhs_t=g%rhs1_t, bcs_t=wrk2d(:, 2))
                case (7)
                    call MatMul_7d_antisym(kmax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), &
                                           u, du1_n, g%periodic, &
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
                        wrk2d(1:len, 2) = wrk2d(1:len, 2) + g%lu1(kmax, idl - ic)*du1_n(:, kmax - ic)
                    end do
                    print *, u(:, kmax), wrk2d(1:len, 2)
                end if

            end do

        end do

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
