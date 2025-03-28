#include "dns_error.h"

module FDM_Derivative
    use TLab_Constants, only: wp, wi, pi_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use FDM_MatMul
    use FDM_PROCS
    use FDM_ComX_Direct
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use FDM_Com0_Jacobian
    implicit none
    private

    type, public :: fdm_derivative_dt
        sequence
        integer mode_fdm                        ! finite-difference method
        integer(wi) size                        ! number of grid points, for convenience in the code
        logical :: periodic = .false.
        logical :: need_1der = .false.          ! In Jacobian formulation, I need 1. order derivative for the 2. order if non-uniform
        integer nb_diag(2)                      ! # of left and right diagonals  (max 5/7)
        real(wp) :: rhs_b(4, 7), rhs_t(4, 7)    ! RHS data for Neumann boundary conditions, max. # of diagonals is 7, # rows is 7/2+1
        real(wp), allocatable :: lhs(:, :)      ! memory space for LHS
        real(wp), allocatable :: rhs(:, :)      ! memory space for RHS
        real(wp), allocatable :: mwn(:)         ! memory space for modified wavenumbers
        real(wp), allocatable :: lu(:, :)       ! memory space for LU decomposition
    end type fdm_derivative_dt

    public :: FDM_Der1_Initialize
    ! public :: FDM_Der1_CreateSystem
    public :: FDM_Der1_Solve

    public :: FDM_Der2_Initialize
    ! public :: FDM_Der2_CreateSystem
    public :: FDM_Der2_Solve

    integer, parameter, public :: FDM_COM4_JACOBIAN = 4
    integer, parameter, public :: FDM_COM6_JACOBIAN_PENTA = 5
    integer, parameter, public :: FDM_COM6_JACOBIAN = 6
    integer, parameter, public :: FDM_COM6_JACOBIAN_HYPER = 7
    integer, parameter, public :: FDM_COM8_JACOBIAN = 8

    integer, parameter, public :: FDM_COM6_DIRECT = 16
    integer, parameter, public :: FDM_COM4_DIRECT = 17

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_Initialize(x, dx, g, periodic, bcs_cases)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:)                   ! Jacobian
        type(fdm_derivative_dt), intent(inout) :: g
        logical, intent(in) :: periodic
        integer, intent(in) :: bcs_cases(:)

        ! -------------------------------------------------------------------
        integer(wi) ib, ip
        integer(wi) nmin, nmax, nsize

        ! ###################################################################
        call FDM_Der1_CreateSystem(x, dx, g, periodic)

        ! -------------------------------------------------------------------
        ! LU decomposition
        if (allocated(g%lu)) deallocate (g%lu)
        if (g%periodic) then
            allocate (g%lu(g%size, g%nb_diag(1) + 2))
        else
            allocate (g%lu(g%size, 5*4))          ! 4 bcs
        end if

        if (periodic) then
            g%lu(:, 1:g%nb_diag(1)) = g%lhs(:, 1:g%nb_diag(1))

            select case (g%nb_diag(1))
            case (3)
                call TRIDPFS(g%size, g%lu(1, 1), g%lu(1, 2), g%lu(1, 3), g%lu(1, 4), g%lu(1, 5))
            case (5)
                call PENTADPFS(g%size, g%lu(1, 1), g%lu(1, 2), g%lu(1, 3), g%lu(1, 4), g%lu(1, 5), g%lu(1, 6), g%lu(1, 7))
            end select

        else                            ! biased,  different BCs
            do ib = 1, size(bcs_cases)
                ip = (ib - 1)*5

                g%lu(:, ip + 1:ip + g%nb_diag(1)) = g%lhs(:, 1:g%nb_diag(1))

                call FDM_Bcs_Neumann(bcs_cases(ib), g%lu(:, ip + 1:ip + g%nb_diag(1)), g%rhs(:, 1:g%nb_diag(2)), g%rhs_b, g%rhs_t)

                nmin = 1; nmax = g%size
                if (any([BCS_ND, BCS_NN] == bcs_cases(ib))) nmin = nmin + 1
                if (any([BCS_DN, BCS_NN] == bcs_cases(ib))) nmax = nmax - 1
                nsize = nmax - nmin + 1

                select case (g%nb_diag(1))
                case (3)
                    call TRIDFS(nsize, g%lu(nmin:, ip + 1), g%lu(nmin:, ip + 2), g%lu(nmin:, ip + 3))
                case (5)
                    call PENTADFS2(nsize, g%lu(nmin:, ip + 1), g%lu(nmin:, ip + 2), g%lu(nmin:, ip + 3), g%lu(nmin:, ip + 4), g%lu(nmin:, ip + 5))
                end select

            end do

        end if

        return
    end subroutine FDM_Der1_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_CreateSystem(x, dx, g, periodic)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(:)                   ! Jacobian
        type(fdm_derivative_dt), intent(inout) :: g
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx
        integer, parameter :: ndl_max = 5, ndr_max = 7

        ! ###################################################################
        g%size = size(x)                ! # grid points
        nx = g%size                     ! for code readability

        if (allocated(g%lhs)) deallocate (g%lhs)
        if (allocated(g%rhs)) deallocate (g%rhs)
        if (allocated(g%mwn)) deallocate (g%mwn)
        allocate (g%lhs(nx, ndl_max))
        allocate (g%rhs(nx, ndr_max))
        allocate (g%mwn(nx))

        g%periodic = periodic

        ! -------------------------------------------------------------------
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)

        case (FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER)
            call FDM_C1N6_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (periodic) then
            nx = g%size

#define wn(i) g%mwn(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn(:) = 2.0_wp*(coef(3)*sin(wn(:)) + coef(4)*sin(2.0_wp*wn(:)) + coef(5)*sin(3.0_wp*wn(:))) &
                            /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(wn(:)))

#undef wn

        end if

        return
    end subroutine FDM_Der1_CreateSystem

! ###################################################################
! ###################################################################
    subroutine FDM_Der1_Solve(nlines, bcs, g, lu1, u, result, wrk2d)
        integer(wi), intent(in) :: nlines   ! # of lines to be solved
        integer(wi), intent(in) :: bcs(2)   ! BCs at xmin (1) and xmax (2):
        !                                   0 biased, non-zero
        !                                   1 forced to zero
        type(fdm_derivative_dt), intent(in) :: g
        real(wp), intent(in) :: lu1(:, :)
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: wrk2d(*)

        integer(wi) nmin, nmax, nsize, ip, ibc

! ###################################################################
        ibc = bcs(1) + bcs(2)*2
        ip = ibc*5

        nmin = 1; nmax = g%size
        if (any([BCS_ND, BCS_NN] == ibc)) then
            result(:, 1) = 0.0_wp      ! homogeneous bcs
            nmin = nmin + 1
        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            result(:, g%size) = 0.0_wp
            nmax = nmax - 1
        end if
        nsize = nmax - nmin + 1

        if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm)) then
            select case (g%nb_diag(2))
            case (3)
                call MatMul_3d(g%rhs(:, 1:3), u, result)
            case (5)
                call MatMul_5d(g%rhs(:, 1:5), u, result)
            end select
        else
            select case (g%nb_diag(2))
            case (3)
                call MatMul_3d_antisym(g%size, nlines, g%rhs(:, 1), g%rhs(:, 2), g%rhs(:, 3), &
                                       u, result, g%periodic, ibc, g%rhs_b, g%rhs_t)
            case (5)
                call MatMul_5d_antisym(g%size, nlines, g%rhs(:, 1), g%rhs(:, 2), g%rhs(:, 3), g%rhs(:, 4), g%rhs(:, 5), &
                                       u, result, g%periodic, ibc, g%rhs_b, g%rhs_t)
            case (7)
                call MatMul_7d_antisym(g%size, nlines, g%rhs(:, 1), g%rhs(:, 2), g%rhs(:, 3), g%rhs(:, 4), g%rhs(:, 5), g%rhs(:, 6), g%rhs(:, 7), &
                                       u, result, g%periodic, ibc, g%rhs_b, g%rhs_t)
            end select
        end if

        if (g%periodic) then
            select case (g%nb_diag(1))
            case (3)
                call TRIDPSS(g%size, nlines, lu1(1, 1), lu1(1, 2), lu1(1, 3), lu1(1, 4), lu1(1, 5), &
                             result, wrk2d)
            case (5)
                call PENTADPSS(g%size, nlines, lu1(1, 1), lu1(1, 2), lu1(1, 3), lu1(1, 4), lu1(1, 5), lu1(1, 6), lu1(1, 7), &
                               result)
            end select

        else
            select case (g%nb_diag(1))
            case (3)
                call TRIDSS(nsize, nlines, lu1(nmin:, ip + 1), lu1(nmin:, ip + 2), lu1(nmin:, ip + 3), &
                            result(:, nmin:))
            case (5)
                call PENTADSS2(nsize, nlines, lu1(nmin:, ip + 1), lu1(nmin:, ip + 2), lu1(nmin:, ip + 3), lu1(nmin:, ip + 4), lu1(nmin:, ip + 5), &
                               result(:, nmin:))
            end select

        end if

        return
    end subroutine FDM_Der1_Solve

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der2_Initialize(x, dx, g, periodic, uniform)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(inout) :: dx(:, :)             ! Jacobians
        type(fdm_derivative_dt), intent(inout) :: g     ! fdm plan for 2. order derivative
        logical, intent(in) :: periodic, uniform

        ! ###################################################################
        call FDM_Der2_CreateSystem(x, dx, g, periodic, uniform)

        ! -------------------------------------------------------------------
        ! LU decomposition
        if (allocated(g%lu)) deallocate (g%lu)
        if (g%periodic) then
            allocate (g%lu(g%size, g%nb_diag(1) + 2))
        else
            allocate (g%lu(g%size, g%nb_diag(1)*1))          ! Only 1 bcs
        end if

        g%lu(:, 1:g%nb_diag(1)) = g%lhs(:, 1:g%nb_diag(1))
        if (g%periodic) then
            select case (g%nb_diag(1))
            case (3)
                call TRIDPFS(g%size, g%lu(1, 1), g%lu(1, 2), g%lu(1, 3), g%lu(1, 4), g%lu(1, 5))
            end select

        else
            select case (g%nb_diag(1))
            case (3)
                call TRIDFS(g%size, g%lu(:, 1), g%lu(:, 2), g%lu(:, 3))
            end select

        end if

        return
    end subroutine FDM_Der2_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der2_CreateSystem(x, dx, g, periodic, uniform)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(inout) :: dx(:, :)             ! Jacobians
        type(fdm_derivative_dt), intent(inout) :: g     ! fdm plan for 2. order derivative
        logical, intent(in) :: periodic, uniform

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx
        integer, parameter :: ndl_max = 5, ndr_max = 7

        ! ###################################################################
        g%size = size(x)                ! # grid points
        nx = g%size                     ! for code readability

        if (allocated(g%lhs)) deallocate (g%lhs)
        if (allocated(g%rhs)) deallocate (g%rhs)
        if (allocated(g%mwn)) deallocate (g%mwn)
        allocate (g%lhs(nx, ndl_max))
        allocate (g%rhs(nx, ndr_max + ndl_max))     ! ndl_max is space for du correction in nonuniform case
        allocate (g%mwn(nx))

        g%periodic = periodic

        ! -------------------------------------------------------------------
        select case (g%mode_fdm)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)
            if (.not. uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
            call FDM_C2N6_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)
            if (.not. uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(g%size, dx, g%lhs, g%rhs, g%nb_diag, coef, periodic)
            if (.not. uniform) g%need_1der = .true.

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)
            g%need_1der = .false.

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(g%size, x, g%lhs, g%rhs, g%nb_diag)
            g%need_1der = .false.

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (periodic) then

#define wn(i) g%mwn(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn(:) = 2.0_wp*(coef(3)*(1.0_wp - cos(wn(:))) + coef(4)*(1.0_wp - cos(2.0_wp*wn(:))) + coef(5)*(1.0_wp - cos(3.0_wp*wn(:)))) &
                       /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(2.0_wp*wn(:)))

#undef wn

        end if

        return
    end subroutine FDM_Der2_CreateSystem

    ! ###################################################################################
    ! ###################################################################################
    subroutine FDM_Der2_Solve(nlines, g, lu, u, result, du, wrk2d)
        integer(wi), intent(in) :: nlines                   ! # of lines to be solved
        type(fdm_derivative_dt), intent(in) :: g            ! plan for 2. order derivative
        real(wp), intent(in) :: lu(:, :)
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(in) :: du(nlines, g%size)          ! 1. derivative for correction in case of Jacobian formulation
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(out) :: wrk2d(*)

        integer(wi) ip

        ! ###################################################################
        if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm)) then
            select case (g%nb_diag(2))
            case (5)
                call MatMul_5d(g%rhs(:, 1:5), u, result)
            end select
        else
            select case (g%nb_diag(2))
            case (5)
                call MatMul_5d_sym(g%size, nlines, g%rhs(:, 1), g%rhs(:, 2), g%rhs(:, 3), g%rhs(:, 4), g%rhs(:, 5), &
                                   u, result, g%periodic)
            case (7)
                call MatMul_7d_sym(g%size, nlines, g%rhs(:, 1), g%rhs(:, 2), g%rhs(:, 3), g%rhs(:, 4), g%rhs(:, 5), g%rhs(:, 6), g%rhs(:, 7), &
                                   u, result, g%periodic)
            end select
        end if

        if (g%need_1der) then
            ip = g%nb_diag(2)      ! add Jacobian correction A_2 dx2 du
            ! so far, only tridiagonal systems
            call MatMul_3d_add(g%size, nlines, g%rhs(:, ip + 1), g%rhs(:, ip + 2), g%rhs(:, ip + 3), du, result)
        end if

        if (g%periodic) then
            select case (g%nb_diag(1))
            case (3)
                call TRIDPSS(g%size, nlines, lu(1, 1), lu(1, 2), lu(1, 3), lu(1, 4), lu(1, 5), &
                             result, wrk2d)
            end select
        else
            select case (g%nb_diag(1))
            case (3)
                call TRIDSS(g%size, nlines, lu(:, 1), lu(:, 2), lu(:, 3), &
                            result)
            end select
        end if

        return
    end subroutine FDM_Der2_Solve

end module FDM_Derivative
