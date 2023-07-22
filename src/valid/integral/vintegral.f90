#include "dns_const.h"

program VINTEGRAL
    use TLAB_CONSTANTS
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLAB_VARS, only: visc, schmidt
    use TLAB_PROCS
    use TLAB_ARRAYS, only: wrk1d, txc, x!, wrk2d, wrk3d
    use FDM_ComX_Direct
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
    integer(wi) bcs_aux(2, 2)
    real(wp) :: lambda, dummy, wk, x_0!, coef(5)
    integer(wi) :: test_type, ibc, ip
    integer(wi) :: nmin, nmax, nsize

    integer, parameter :: i1 = 1, cases(2) = [1, 2]

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

! Valid stettings
    test_type = 0

    g%periodic = .false.
    wk = 1 ! WRITE(*,*) 'Wavenumber ?'; READ(*,*) wk
    lambda = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda
    g%mode_fdm1 = FDM_COM6_JACOBIAN ! FDM_COM6_JACOBIAN_PENTA
    if (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) C1N6M_ALPHA = 0.56
    g%mode_fdm2 = g%mode_fdm1

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
        u(:, i) = sin(2.0_wp*pi_wp/g%scale*wk*g%nodes(i)) ! + pi_wp/4.0_wp)
        du1_a(:, i) = (2.0_wp*pi_wp/g%scale*wk) &
                      *cos(2.0_wp*pi_wp/g%scale*wk*g%nodes(i))! + pi_wp/4.0_wp)
! Gaussian
        ! u(i)     = EXP(-(g%nodes(i)-x_0*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        ! du1_a(i) =-(g%nodes(i)-x_0*g%scale)/(g%scale/wk)**2*u(i)
        ! v(i)     = EXP(-(g%nodes(i)-x_0*C_05_R*g%scale)**2/(2.0_wp*(g%scale/wk)**2))
        ! dv1_a(i) =-(g%nodes(i)-x_0*C_05_R*g%scale)/(g%scale/wk)**2*v(i)
! exponential
        ! u(i) = EXP(-g%nodes(i)*lambda)
        ! du1_a(i) = -lambda*u(i)
        ! du2_a(i) =  lambda*lambda*u(i)
        ! v(i) = EXP(g%nodes(i)*x_0/g%scale)
        ! dv1_a(i) = x_0/g%scale*v(i)
! step
        ! u(i) = MAX(0.0_wp,(g%nodes(i)-g%nodes(imax/2))*x_0)
        ! du1_a(i) = (1.0_wp+SIGN(1.0_wp,g%nodes(i)-g%nodes(imax/2)))*C_05_R*x_0
! tanh
        ! u(i) = x_0 * LOG( 1.0_wp + EXP( (g%nodes(i)-g%nodes(imax/2))/x_0 ) )
        ! du1_a(i) = C_05_R*( 1.0_wp + TANH( C_05_R*(g%nodes(i)-g%nodes(imax/2))/x_0 ) )
! Polynomial
        ! dummy = 4.0_wp
        ! u(l, i) = ((g%scale - g%nodes(i))/lambda)**dummy
        ! du1_a(l, i) = -dummy*((g%scale - g%nodes(i))/lambda)**(dummy - 1.0_wp)
! zero
        ! u(i) = 0.0_wp
        ! du1_a(i) = 0.0_wp
    end do

! ###################################################################
! Integral
! ###################################################################
    if (test_type == 0) then

        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g, u, f)
        f = du1_a

        do ip = 1, size(cases)
            ibc = cases(ip)
            print *, new_line('a'), 'Case ', ibc

            nmin = 1
            nmax = imax
            select case (ibc)
            case (1)
                nmin = nmin + 1
                bcs(:, 1) = u(:, 1)
            case (2)
                nmax = nmax - 1
                bcs(:, 1) = u(:, imax)
            end select
            nsize = nmax - nmin + 1

            if (g%mode_fdm1 == FDM_COM6_JACOBIAN) then
                call INT_C1N6_LHS(imax, ibc, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
                call INT_C1N6_RHS(imax, len, ibc, g%jac, f, w_n)

            elseif (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
                call INT_C1N6M_LHS(imax, ibc, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                call INT_C1N6M_RHS(imax, len, ibc, g%jac, f, w_n)

            end if

            if (g%mode_fdm1 == FDM_COM6_JACOBIAN) then
                call PENTADFS(nsize, wrk1d(nmin:nmax, 1), wrk1d(nmin:nmax, 2), wrk1d(nmin:nmax, 3), wrk1d(nmin:nmax, 4), wrk1d(nmin:nmax, 5))
 call PENTADSS(nsize, len, wrk1d(nmin:nmax, 1), wrk1d(nmin:nmax, 2), wrk1d(nmin:nmax, 3), wrk1d(nmin:nmax, 4), wrk1d(nmin:nmax, 5), w_n(:, nmin:nmax))

            elseif (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
                call HEPTADFS(nsize, wrk1d(nmin:nmax, 1), wrk1d(nmin:nmax, 2), wrk1d(nmin:nmax, 3), wrk1d(nmin:nmax, 4), wrk1d(nmin:nmax, 5), wrk1d(nmin:nmax, 6), wrk1d(nmin:nmax, 7))
                call HEPTADSS(nsize, len, wrk1d(nmin:nmax, 1), wrk1d(nmin:nmax, 2), wrk1d(nmin:nmax, 3), wrk1d(nmin:nmax, 4), wrk1d(nmin:nmax, 5), wrk1d(nmin:nmax, 6), wrk1d(nmin:nmax, 7), w_n(:, nmin:nmax))

            end if

            select case (ibc)
            case (1)                    ! BCs at the bottom
                w_n(:, 1) = 0.0_wp
            case (2)                    ! BCs at the top
                w_n(:, imax) = 0.0_wp
            end select
            do i = 1, imax
                w_n(:, i) = w_n(:, i) + bcs(:, 1)
            end do

            call check(u, w_n, 'integral.dat')

        end do

! ! ###################################################################
! ! First order equation
! ! ###################################################################
!     else if (test_type == 1) then

!         f = du1_a + lambda*u

!         ibc = 2

!         if (g%mode_fdm1 == FDM_COM6_JACOBIAN) then
!             call INT_C1N6_LHS_E(imax, ibc, g%jac, lambda, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
!             call INT_C1N6_RHS(imax, i1, ibc, g%jac, f, w_n)
!         elseif (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
! call INT_C1N6M_LHS_E(imax, ibc, g%jac, lambda, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7), wrk1d(1, 8))
!             call INT_C1N6M_RHS(imax, i1, ibc, g%jac, f, w_n)
!         end if

!         if (ibc == 1) then
!             if (g%mode_fdm1 == FDM_COM6_JACOBIAN) then
!                 call PENTADFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
!                 call PENTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), w_n(2))
!                 call PENTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
!                 dummy = w_n(1); w_n(1) = 0.0_wp
!                 w_n = w_n + u(1)*wrk1d(1:imax, 6) ! BCs
!                 dummy = (dummy + wrk1d(1, 3)*w_n(1) + wrk1d(1, 4)*w_n(2) + wrk1d(1, 5)*w_n(3))/g%jac(1, 1)
!             elseif (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!                 call HEPTADFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7))
!                 call HEPTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7), w_n(2))
!                 call HEPTADSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6), wrk1d(2, 7), wrk1d(2, 8))
!                 dummy = w_n(1); w_n(1) = 0.0_wp
!                 w_n = w_n + u(1)*wrk1d(1:imax, 8) ! BCs
!                 dummy = (dummy + wrk1d(1, 4)*w_n(1) + wrk1d(1, 5)*w_n(2) + wrk1d(1, 6)*w_n(3))/g%jac(1, 1)
!             end if

!         else if (ibc == 2) then
!             if (g%mode_fdm1 == FDM_COM6_JACOBIAN) then
!                 call PENTADFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
!                 call PENTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), w_n(1))
!                 call PENTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
!                 dummy = w_n(imax); w_n(imax) = 0.0_wp
!                 w_n = w_n + u(imax)*wrk1d(1:imax, 6) ! BCs
!                 dummy = (dummy + wrk1d(imax, 1)*w_n(imax - 2) + wrk1d(imax, 2)*w_n(imax - 1) + wrk1d(imax, 3)*w_n(imax))/g%jac(imax, 1)
!             elseif (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!                 call HEPTADFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
!                 call HEPTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7), w_n(1))
!                 call HEPTADSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7), wrk1d(1, 8))
!                 dummy = w_n(imax); w_n(imax) = 0.0_wp
!                 w_n = w_n + u(imax)*wrk1d(1:imax, 8) ! BCs
!                 dummy = (dummy + wrk1d(imax, 2)*w_n(imax - 2) + wrk1d(imax, 3)*w_n(imax - 1) + wrk1d(imax, 4)*w_n(imax))/g%jac(imax, 1)
!             end if
!         end if

!         if (ibc == 1) then
!             write (*, *) dummy, dw1_n(1)
!         else
!             write (*, *) dummy, dw1_n(imax)
!         end if

!         f = u

! ! ###################################################################
! ! Second order equation
! ! ###################################################################
!     else if (test_type == 2) then

!         call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g, u, f, tmp)

!         dummy = lambda*lambda
!         f = du2_a - dummy*u

!         ! to be rewritten in terms of INT_C2NX_
!         ! call INT_C2N6_LHS_E(imax, g%jac, ibc, dummy, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
!         call PENTADFS(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
!         call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
!         call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
!         ! call INT_C2N6_RHS(imax, i1, g%jac, f, w_n)
!         call PENTADSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), w_n(2))
!         w_n(:) = w_n(:) + u(1)*wrk1d(:, 6) + u(imax)*wrk1d(:, 7)

    end if

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

! ! ###################################################################
! ! IO - Error and function values
!     open (20, file='integral.dat')
!     error = 0.0_wp
!     sol = 0.0_wp
!     imin_loc = 1; imax_loc = imax
!     do i = imin_loc, imax_loc
!         write (20, 1000) g%nodes(i), u(i), w_n(i), u(i) - w_n(i)
!         w_n(i) = abs(u(i) - w_n(i))
!         error = error + w_n(i)*w_n(i)
!         sol = sol + u(i)*u(i)
!     end do
!     close (20)

!     write (*, *) 'Solution L2-norm ......:', sqrt(g%jac(1, 1)*sol)
!     write (*, *) 'Error L2-norm .........:', sqrt(g%jac(1, 1)*error)
!     write (*, *) 'Error Linf-norm .......:', maxval(w_n(1:imax))
!     write (*, *) 'Relative error ........:', sqrt(error)/sqrt(sol)
!     write (*, *) 'Derivative overshoot ..:', minval(dw1_n(1:imax))

!     stop

! 1000 format(6(1x, e17.10e3))

end program VINTEGRAL
