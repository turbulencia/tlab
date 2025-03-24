#include "dns_error.h"

module FDM_Derivative
    use TLab_Constants, only: wp, wi, pi_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use FDM_MatMul
    use FDM_ComX_Direct
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use FDM_Com0_Jacobian
    use FDM_Interpolate, only: fdm_interpol_dt
    implicit none
    private

    ! This type should be split into 2 derivative types and the interpolation data moved to only fdm.
    type, public :: fdm_dt
        sequence
        character*8 name
        integer(wi) size
        logical :: uniform = .false.
        logical :: periodic = .false.
        real(wp) scale
        real(wp), pointer :: memory(:, :)       ! memory space
        !
        real(wp), pointer :: nodes(:)
        real(wp), pointer :: jac(:, :)          ! pointer to Jacobians
        !
        integer mode_fdm1                       ! finite-difference method for 1. order derivative
        integer nb_diag_1(2)                    ! # of left and right diagonals 1. order derivative (max 5/7)
        real(wp) :: rhs1_b(4, 7), rhs1_t(4, 7)  ! RHS data for Neumann boundary conditions, 1. order derivative max. # of diagonals is 7, # rows is 7/2+1
        real(wp), pointer :: lhs1(:, :)         ! pointer to LHS for 1. derivative
        real(wp), pointer :: rhs1(:, :)         ! pointer to RHS for 1. derivative
        real(wp), pointer :: mwn1(:)            ! pointer to modified wavenumbers
        real(wp), pointer :: lu1(:, :)          ! pointer to LU decomposition for 1. derivative
        !
        integer mode_fdm2                       ! finite-difference method for 2. order derivative
        integer nb_diag_2(2)                    ! # of left and right diagonals 2. order derivative (max 5/7)
        real(wp), pointer :: lhs2(:, :)         ! pointer to LHS for 2. derivative
        real(wp), pointer :: rhs2(:, :)         ! pointer to RHS for 2. derivative
        real(wp), pointer :: mwn2(:)            ! pointer to modified wavenumbers
        real(wp), pointer :: lu2(:, :)          ! pointer to LU decomposition for 2. derivative
        logical :: need_1der = .false.          ! In Jacobian formulation, I need 1. order derivative for the 2. order if non-uniform
        !
        type(fdm_interpol_dt) :: intl

    end type fdm_dt

    public :: FDM_Der1_CreateSystem
    public :: FDM_Der1_Solve

    public :: FDM_Der2_CreateSystem
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
    subroutine FDM_Der1_CreateSystem(g, periodic)
        type(fdm_dt), intent(inout) :: g
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx

        ! ###################################################################
        select case (g%mode_fdm1)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(g%size, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, periodic)

        case (FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER)
            call FDM_C1N6_Jacobian(g%size, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, periodic)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(g%size, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(g%size, g%nodes, g%lhs1, g%rhs1, g%nb_diag_1)

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(g%size, g%nodes, g%lhs1, g%rhs1, g%nb_diag_1)

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (periodic) then
            nx = g%size

#define wn(i) g%mwn1(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn1(:) = 2.0_wp*(coef(3)*sin(wn(:)) + coef(4)*sin(2.0_wp*wn(:)) + coef(5)*sin(3.0_wp*wn(:))) &
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
        type(fdm_dt), intent(in) :: g
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

        if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm1)) then
            select case (g%nb_diag_1(2))
            case (3)
                call MatMul_3d(g%rhs1(:, 1:3), u, result)
            case (5)
                call MatMul_5d(g%rhs1(:, 1:5), u, result)
            end select
        else
            select case (g%nb_diag_1(2))
            case (3)
                call MatMul_3d_antisym(g%size, nlines, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), &
                                       u, result, g%periodic, ibc, g%rhs1_b, g%rhs1_t)
            case (5)
                call MatMul_5d_antisym(g%size, nlines, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), &
                                       u, result, g%periodic, ibc, g%rhs1_b, g%rhs1_t)
            case (7)
                call MatMul_7d_antisym(g%size, nlines, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), &
                                       u, result, g%periodic, ibc, g%rhs1_b, g%rhs1_t)
            end select
        end if

        if (g%periodic) then
            select case (g%nb_diag_1(1))
            case (3)
                call TRIDPSS(g%size, nlines, lu1(1, 1), lu1(1, 2), lu1(1, 3), lu1(1, 4), lu1(1, 5), &
                             result, wrk2d)
            case (5)
                call PENTADPSS(g%size, nlines, lu1(1, 1), lu1(1, 2), lu1(1, 3), lu1(1, 4), lu1(1, 5), lu1(1, 6), lu1(1, 7), &
                               result)
            end select

        else
            select case (g%nb_diag_1(1))
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
    subroutine FDM_Der2_CreateSystem(g, periodic)
        type(fdm_dt), intent(inout) :: g
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx

        ! ###################################################################
        select case (g%mode_fdm2)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(g%size, g%jac(:, 2), g%lhs2, g%rhs2, g%nb_diag_2, coef, periodic)
            if (.not. g%uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
            call FDM_C2N6_Jacobian(g%size, g%jac(:, 2), g%lhs2, g%rhs2, g%nb_diag_2, coef, periodic)
            if (.not. g%uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(g%size, g%jac(:, 2), g%lhs2, g%rhs2, g%nb_diag_2, coef, periodic)
            if (.not. g%uniform) g%need_1der = .true.

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(g%size, g%nodes, g%lhs2, g%rhs2, g%nb_diag_2)

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(g%size, g%nodes, g%lhs2, g%rhs2, g%nb_diag_2)

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (periodic) then
            nx = g%size

#define wn(i) g%mwn2(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn2(:) = 2.0_wp*(coef(3)*(1.0_wp - cos(wn(:))) + coef(4)*(1.0_wp - cos(2.0_wp*wn(:))) + coef(5)*(1.0_wp - cos(3.0_wp*wn(:)))) &
                        /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(2.0_wp*wn(:)))

#undef wn

        end if

        return
    end subroutine FDM_Der2_CreateSystem

    ! ###################################################################################
    ! ###################################################################################
    subroutine FDM_Der2_Solve(nlines, g, lu2, u, result, du, wrk2d)
        integer(wi), intent(in) :: nlines                  ! # of lines to be solved
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: lu2(:, :)
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(in) :: du(nlines, g%size)          ! 1. derivative for correction in case of Jacobian formulation
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(out) :: wrk2d(*)

        integer(wi) ip

        ! ###################################################################
        if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm2)) then
            select case (g%nb_diag_2(2))
            case (5)
                call MatMul_5d(g%rhs2(:, 1:5), u, result)
            end select
        else
            select case (g%nb_diag_2(2))
            case (5)
                call MatMul_5d_sym(g%size, nlines, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), &
                                   u, result, g%periodic)
            case (7)
                call MatMul_7d_sym(g%size, nlines, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), g%rhs2(:, 6), g%rhs2(:, 7), &
                                   u, result, g%periodic)
            end select
            if (g%need_1der) then
                ip = g%nb_diag_2(2)      ! add Jacobian correction A_2 dx2 du
                ! so far, only tridiagonal systems
                call MatMul_3d_add(g%size, nlines, g%rhs2(:, ip + 1), g%rhs2(:, ip + 2), g%rhs2(:, ip + 3), du, result)
            end if
        end if

        if (g%periodic) then
            select case (g%nb_diag_2(1))
            case (3)
                call TRIDPSS(g%size, nlines, lu2(1, 1), lu2(1, 2), lu2(1, 3), lu2(1, 4), lu2(1, 5), &
                             result, wrk2d)
            end select
        else
            select case (g%nb_diag_2(1))
            case (3)
                call TRIDSS(g%size, nlines, lu2(:, 1), lu2(:, 2), lu2(:, 3), &
                            result)
            end select
        end if

        return
    end subroutine FDM_Der2_Solve

end module FDM_Derivative
