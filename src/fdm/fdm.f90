#include "dns_error.h"

module FDM
    use TLab_Constants, only: wp, wi
    use FDM_Integral
    implicit none
    private

    type, public :: fdm_dt
        sequence
        character*8 name
        integer(wi) size
        logical uniform, periodic
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
        real(wp), pointer :: lu1(:, :)          ! pointer to LU decomposition for 1. derivative
        real(wp), pointer :: mwn1(:)            ! pointer to modified wavenumbers
        !
        real(wp) :: rhsr_b(5, 0:7), rhsr_t(0:4, 8)      ! RHS data for reduced boundary conditions; max. # of diagonals is 7, # rows is 7/2+1
        !
        type(fdm_integral_dt) :: fdmi(2)
        !
        real(wp), pointer :: lu0i(:, :)                 ! pointer to LU decomposition for interpolation
        real(wp), pointer :: lu1i(:, :)                 ! pointer to LU decomposition for 1. derivative inc. interp.
        !
        integer mode_fdm2                   ! finite-difference method for 2. order derivative
        integer nb_diag_2(2)                ! # of left and right diagonals 2. order derivative (max 5/7)
        real(wp), pointer :: lhs2(:, :)     ! pointer to LHS for 2. derivative
        real(wp), pointer :: rhs2(:, :)     ! pointer to RHS for 2. derivative
        real(wp), pointer :: lu2(:, :)      ! pointer to LU decomposition for 2. derivative
        real(wp), pointer :: mwn2(:)        ! pointer to modified wavenumbers
        logical :: need_1der = .false.      ! In Jacobian formulation, I need 1. order derivative for the 2. order if non-uniform

    end type fdm_dt

    type(fdm_dt), public :: g(3)                ! Grid information along 3 directions

    public :: FDM_Initialize

    integer, parameter, public :: FDM_COM4_JACOBIAN = 4
    integer, parameter, public :: FDM_COM6_JACOBIAN_PENTA = 5
    integer, parameter, public :: FDM_COM6_JACOBIAN = 6
    integer, parameter, public :: FDM_COM6_JACOBIAN_HYPER = 7
    integer, parameter, public :: FDM_COM8_JACOBIAN = 8

    integer, parameter, public :: FDM_COM6_DIRECT = 16
    integer, parameter, public :: FDM_COM4_DIRECT = 17

contains
    subroutine FDM_Initialize(g, nodes, locScale)
        use TLab_Constants, only: pi_wp, efile, wfile, BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, roundoff_wp
#ifdef TRACE_ON
        use TLab_Constants, only: tfile
#endif
        use TLab_WorkFlow, only: stagger_on
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        use TLab_Memory, only: TLab_Allocate_Real
        use FDM_PROCS, only: FDM_Bcs_Neumann
        use FDM_MatMul
        use FDM_ComX_Direct
        use FDM_Com1_Jacobian
        use FDM_Com2_Jacobian

        type(fdm_dt), intent(inout) :: g                    ! grid structure
        real(wp), intent(in) :: nodes(:)                    ! positions of the grid nodes
        real(wp), intent(in), optional :: locScale          ! for consistency check

! -------------------------------------------------------------------
        integer(wi) i, ib, ip, ig, nx, ndl, ndr, inb_grid
        integer(wi) nmin, nmax, nsize, bcs_cases(4)
        real(wp) coef(5)

        integer, parameter :: i1 = 1

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'Entering '//__FILE__)
#endif

        ! ###################################################################
        ! Consistency check
        ! ###################################################################
        if (g%periodic .and. g%mode_fdm1 == FDM_COM4_DIRECT) g%mode_fdm1 = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm1 == FDM_COM6_DIRECT) g%mode_fdm1 = FDM_COM6_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm2 == FDM_COM4_DIRECT) g%mode_fdm2 = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm2 == FDM_COM6_DIRECT) g%mode_fdm2 = FDM_COM6_JACOBIAN        ! they are the same for uniform grids.

        if (any([FDM_COM4_DIRECT, FDM_COM6_DIRECT] == g%mode_fdm1)) g%mode_fdm1 = FDM_COM6_JACOBIAN ! undeveloped; I would need to read separately 1. and 2. order information
        if (any([FDM_COM6_JACOBIAN_PENTA] == g%mode_fdm2)) g%mode_fdm2 = FDM_COM6_JACOBIAN          ! undeveloped; I would need to read separately 1. and 2. order information
        if (g%mode_fdm2 == FDM_COM6_JACOBIAN) g%mode_fdm2 = FDM_COM6_JACOBIAN_HYPER                 ! default

        if (g%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then                                            ! CFL_max depends on max[g%mwn1(:)]
            call TLab_Write_ASCII(wfile, __FILE__//'. Main.SpaceOrder.CompactJacobian6Penta requires adjusted CFL-number depending on alpha and beta values.')
        end if

        if (size(nodes) /= g%size) then
            call TLab_Write_ASCII(efile, __FILE__//'. Unmatched grid size.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (g%size > 1) then
            g%scale = nodes(g%size) - nodes(1)
            if (g%periodic) g%scale = g%scale*(1.0_wp + 1.0_wp/real(g%size - 1, wp))
        else
            g%scale = 1.0_wp  ! to avoid conditionals and NaN in some of the calculations below
        end if

        if (present(locScale)) then
! print *, abs((scale_loc - g%scale)/scale_loc)
            if (abs((locScale - g%scale)/g%scale) > roundoff_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Unmatched domain scale.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        ! ###################################################################
        ! Memory allocation
        ! ###################################################################
        inb_grid = 1                            ! Nodes
        inb_grid = inb_grid &
                   + 2 &                        ! Jacobians of first- and second-order derivatives
                   + 2                          ! 1/dx and 1/dx**2 used in time-step stability constraint

        inb_grid = inb_grid &
                   + 5 &                        ! max # of diagonals in LHS for 1. order derivative
                   + 7                          ! max # of diagonals in RHS for 1. order derivative
        if (g%periodic) then
            inb_grid = inb_grid &
                       + 5 + 2 &                ! LU decomposition 1. order
                       + 1                      ! modified wavenumbers for 1. order derivative
        else
            inb_grid = inb_grid &
                       + 5*4                    ! LU decomposition 1. order, 4 bcs
        end if

        inb_grid = inb_grid &
                   + 5 &                        ! max # of diagonals in LHS for 2. order derivative
                   + 7 + 5                      ! max # of diagonals in RHS for 2. order + diagonals for Jacobian case
        if (g%periodic) then
            inb_grid = inb_grid &
                       + 5 + 2 &                ! LU decomposition 2. order
                       + 1                      ! modified wavenumbers for 2. order derivative
        else
            inb_grid = inb_grid &
                       + 5                      ! LU decomposition 2. order, 1 bcs
        end if

        ! inb_grid = inb_grid &
        !            + 5*2 &                      ! max # of diagonals in LHS for 1. integral, 2 bcs
        !            + 7*2                        ! max # of diagonals in RHS for 1. integral, 2 bcs

        if (stagger_on .and. g%periodic) then
            inb_grid = inb_grid &
                       + 5 &                    ! LU decomposition interpolation
                       + 5                      ! LU decomposition 1. order with interpolation
        end if

        ! call TLab_Allocate_Real(__FILE__, g%memory, [g%size, inb_grid], g%name)
        allocate (g%memory(g%size, inb_grid))

        ! ###################################################################
        ! Setting pointers and filling FDM data
        ! ###################################################################
        nx = g%size                     ! node number, for clarity below

        ig = 1                          ! Initialize counter to define pointers inside array x

        ! ###################################################################
        ! Node positions
        ! ###################################################################
        g%nodes => g%memory(:, ig)             ! Define pointer inside x

        g%nodes(:) = nodes(1:nx)        ! Calculate data

        ig = ig + 1                     ! Advance counter

        ! ###################################################################
        ! Space for Jacobians: computational grid is uniform
        ! ###################################################################
        g%jac => g%memory(:, ig:)

        if (nx == 1) then
            g%jac(:, :) = 1.0_wp
            return
        end if

        ig = ig + 4

        ! ###################################################################
        ! first-order derivative
        ! ###################################################################
        g%lhs1 => g%memory(:, ig:)
        ig = ig + 5
        g%rhs1 => g%memory(:, ig:)
        ig = ig + 7

        ! -------------------------------------------------------------------
        ! Jacobian; computational grid is uniform; used for the stencils below and also as grid spacing in the code
        g%jac(:, 1) = 1.0_wp

        select case (g%mode_fdm1)
        case (FDM_COM4_JACOBIAN, FDM_COM4_DIRECT)
            call FDM_C1N4_Jacobian(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef)
            call MatMul_3d_antisym(nx, 1, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%nodes(:), g%jac(:, 1), periodic=.false.)

        case (FDM_COM6_JACOBIAN, FDM_COM6_DIRECT)
            call FDM_C1N6_Jacobian(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef)
            call MatMul_5d_antisym(nx, 1, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%nodes(:), g%jac(:, 1), periodic=.false.)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef)
            call MatMul_7d_antisym(nx, 1, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), g%rhs1(:, 4), g%rhs1(:, 5), g%rhs1(:, 6), g%rhs1(:, 7), g%nodes(:), g%jac(:, 1), periodic=.false.)

        end select

        select case (g%nb_diag_1(1))
        case (3)
            call TRIDFS(nx, g%lhs1(:, 1), g%lhs1(:, 2), g%lhs1(:, 3))
            call TRIDSS(nx, i1, g%lhs1(:, 1), g%lhs1(:, 2), g%lhs1(:, 3), g%jac(1, 1))
        case (5)
            call PENTADFS2(nx, g%lhs1(:, 1), g%lhs1(:, 2), g%lhs1(:, 3), g%lhs1(:, 4), g%lhs1(:, 5))
            call PENTADSS2(nx, i1, g%lhs1(:, 1), g%lhs1(:, 2), g%lhs1(:, 3), g%lhs1(:, 4), g%lhs1(:, 5), g%jac(1, 1))
        end select

        ! Saving operations for the time-stability constraint
        g%jac(:, 3) = 1.0_wp/g%jac(:, 1)
        g%jac(:, 4) = g%jac(:, 3)*g%jac(:, 3)

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        select case (g%mode_fdm1)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_Jacobian(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(nx, g%jac, g%lhs1, g%rhs1, g%nb_diag_1, coef, g%periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(nx, g%nodes, g%lhs1, g%rhs1, g%nb_diag_1)

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(nx, g%nodes, g%lhs1, g%rhs1, g%nb_diag_1)

        end select
        ndl = g%nb_diag_1(1)    ! for readability of the source code
        ndr = g%nb_diag_1(2)

        ! -------------------------------------------------------------------
        ! LU decomposition
        g%lu1 => g%memory(:, ig:)

        g%lu1(:, 1:g%nb_diag_1(1)) = g%lhs1(:, 1:g%nb_diag_1(1))
        if (g%periodic) then
            select case (g%nb_diag_1(1))
            case (3)
                call TRIDPFS(nx, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5))
            case (5)
                call PENTADPFS(nx, g%lu1(1, 1), g%lu1(1, 2), g%lu1(1, 3), g%lu1(1, 4), g%lu1(1, 5), g%lu1(1, 6), g%lu1(1, 7))
            end select

            ig = ig + g%nb_diag_1(1) + 2

        else                            ! biased,  different BCs
            bcs_cases(1:4) = [BCS_DD, BCS_ND, BCS_DN, BCS_NN]
            do ib = 1, 4
                ip = (ib - 1)*5

                g%lu1(:, ip + 1:ip + g%nb_diag_1(1)) = g%lhs1(:, 1:g%nb_diag_1(1))

                call FDM_Bcs_Neumann(bcs_cases(ib), g%lu1(:, ip + 1:ip + g%nb_diag_1(1)), g%rhs1(:, 1:g%nb_diag_1(2)), g%rhs1_b, g%rhs1_t)

                nmin = 1; nmax = nx
                if (any([BCS_ND, BCS_NN] == bcs_cases(ib))) nmin = nmin + 1
                if (any([BCS_DN, BCS_NN] == bcs_cases(ib))) nmax = nmax - 1
                nsize = nmax - nmin + 1

                select case (g%nb_diag_1(1))
                case (3)
                    call TRIDFS(nsize, g%lu1(nmin:, ip + 1), g%lu1(nmin:, ip + 2), g%lu1(nmin:, ip + 3))
                case (5)
                    call PENTADFS2(nsize, g%lu1(nmin:, ip + 1), g%lu1(nmin:, ip + 2), g%lu1(nmin:, ip + 3), g%lu1(nmin:, ip + 4), g%lu1(nmin:, ip + 5))
                end select

                ig = ig + 5

            end do

        end if

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (g%periodic) then
            g%mwn1 => g%memory(:, ig)

#define wn(i) g%mwn1(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            if (.not. stagger_on) then

                g%mwn1(:) = 2.0_wp*(coef(3)*sin(wn(:)) + coef(4)*sin(2.0_wp*wn(:)) + coef(5)*sin(3.0_wp*wn(:))) &
                            /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(wn(:)))

            else ! staggered case has different modified wavenumbers!

                select case (g%mode_fdm1)

                case DEFAULT
                    coef = [9.0_wp/62.0_wp, 0.0_wp, 63.0_wp/62.0_wp, 17.0_wp/62.0_wp, 0.0_wp]

                end select

                g%mwn1(:) = 2.0_wp*(coef(3)*sin(1.0_wp/2.0_wp*wn(:)) + coef(4)/3.0_wp*sin(3.0_wp/2.0_wp*wn(:))) &
                            /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)))

            end if

#undef wn
            g%mwn1(:) = (g%mwn1(:)/g%jac(1, 1))**2      ! as used in Poisson solver

            ig = ig + 1

        end if

! ###################################################################
! first-order integrals (cases lambda = 0.0_wp)
! ###################################################################
        if (.not. g%periodic) then
            bcs_cases(1:2) = [BCS_MIN, BCS_MAX]
            do ib = 1, 2
                g%fdmi(ib)%bc = bcs_cases(ib)
                call FDM_Int1_Initialize(g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), 0.0_wp, g%fdmi(ib))

                ! LU decomposition
                select case (ndr)
                case (3)
                    call TRIDFS(g%size - 2, g%fdmi(ib)%lhs(2:, 1), g%fdmi(ib)%lhs(2:, 2), g%fdmi(ib)%lhs(2:, 3))
                case (5)
                    call PENTADFS(g%size - 2, g%fdmi(ib)%lhs(2:, 1), g%fdmi(ib)%lhs(2:, 2), g%fdmi(ib)%lhs(2:, 3), &
                                  g%fdmi(ib)%lhs(2:, 4), g%fdmi(ib)%lhs(2:, 5))
                case (7)
                    call HEPTADFS(g%size - 2, g%fdmi(ib)%lhs(2:, 1), g%fdmi(ib)%lhs(2:, 2), g%fdmi(ib)%lhs(2:, 3), &
                                  g%fdmi(ib)%lhs(2:, 4), g%fdmi(ib)%lhs(2:, 5), g%fdmi(ib)%lhs(2:, 6), g%fdmi(ib)%lhs(2:, 7))
                end select
            end do
        end if

! ###################################################################
! second-order derivative: LU factorization done in routine TRID*FS
! ###################################################################
        g%lhs2 => g%memory(:, ig:)
        ig = ig + 5
        g%rhs2 => g%memory(:, ig:)
        ig = ig + 7 + 5

        ! -------------------------------------------------------------------
        ! Jacobian; computational grid is uniform; only used to calculate the stencils in the section below
        g%jac(:, 2) = 1.0_wp

        select case (g%mode_fdm2)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(nx, g%jac(:, 2), g%lhs2, g%rhs2, g%nb_diag_2, coef)
            call MatMul_5d_sym(nx, 1, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), g%nodes(:), g%jac(:, 2), periodic=.false.)

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(nx, g%jac(:, 2), g%lhs2, g%rhs2, g%nb_diag_2, coef)
            call MatMul_5d_sym(nx, 1, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), g%nodes(:), g%jac(:, 2), periodic=.false.)

        case (FDM_COM6_JACOBIAN_HYPER, FDM_COM6_DIRECT, FDM_COM6_JACOBIAN_PENTA)
            call FDM_C2N6_Hyper_Jacobian(nx, g%jac(:, 2), g%lhs2, g%rhs2, g%nb_diag_2, coef)
         call MatMul_7d_sym(nx, 1, g%rhs2(:, 1), g%rhs2(:, 2), g%rhs2(:, 3), g%rhs2(:, 4), g%rhs2(:, 5), g%rhs2(:, 6), g%rhs2(:, 7), g%nodes(:), g%jac(:, 2), periodic=.false.)

        end select

        select case (g%nb_diag_2(1))
        case (3)
            call TRIDFS(nx, g%lhs2(:, 1), g%lhs2(:, 2), g%lhs2(:, 3))
            call TRIDSS(nx, i1, g%lhs2(:, 1), g%lhs2(:, 2), g%lhs2(:, 3), g%jac(1, 2))
        case (5)
            call PENTADFS2(nx, g%lhs2(:, 1), g%lhs2(:, 2), g%lhs2(:, 3), g%lhs2(:, 4), g%lhs2(:, 5))
            call PENTADSS2(nx, i1, g%lhs2(:, 1), g%lhs2(:, 2), g%lhs2(:, 3), g%lhs2(:, 4), g%lhs2(:, 5), g%jac(1, 2))
        end select

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        select case (g%mode_fdm2)

        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(nx, g%jac, g%lhs2, g%rhs2, g%nb_diag_2, coef, g%periodic)
            if (.not. g%uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(nx, g%jac, g%lhs2, g%rhs2, g%nb_diag_2, coef, g%periodic)
            if (.not. g%uniform) g%need_1der = .true.

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(nx, g%jac, g%lhs2, g%rhs2, g%nb_diag_2, coef, g%periodic)
            if (.not. g%uniform) g%need_1der = .true.

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(nx, g%nodes, g%lhs2, g%rhs2, g%nb_diag_2)

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(nx, g%nodes, g%lhs2, g%rhs2, g%nb_diag_2)

        end select

        ! -------------------------------------------------------------------
        ! LU decomposition and wave numbers
        g%lu2 => g%memory(:, ig:)

        g%lu2(:, 1:g%nb_diag_2(1)) = g%lhs2(:, 1:g%nb_diag_2(1))
        if (g%periodic) then
            select case (g%nb_diag_2(1))
            case (3)
                call TRIDPFS(nx, g%lu2(1, 1), g%lu2(1, 2), g%lu2(1, 3), g%lu2(1, 4), g%lu2(1, 5))
            end select
            ig = ig + g%nb_diag_2(1) + 2

        else
            select case (g%nb_diag_2(1))
            case (3)
                call TRIDFS(nx, g%lu2(:, 1), g%lu2(:, 2), g%lu2(:, 3))
            end select
            ig = ig + g%nb_diag_2(1)

        end if

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        if (g%periodic) then
            g%mwn2 => g%memory(:, ig)

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

            g%mwn2(:) = g%mwn2(:)/(g%jac(1, 1)**2)  ! as used in the Helmholtz solver

            ig = ig + 1

        end if

! ###################################################################
! LU factorization interpolation, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
        if ((stagger_on) .and. g%periodic) then
            g%lu0i => g%memory(:, ig:)

            select case (g%mode_fdm1)
            case DEFAULT
                call FDM_C0INT6P_LHS(nx, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3))
            end select
            call TRIDPFS(nx, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5))
            ig = ig + 5
        end if

! ###################################################################
! LU factorization first interp. derivative, done in routine TRID*FS
! ###################################################################
! -------------------------------------------------------------------
! Periodic case; pentadiagonal
! -------------------------------------------------------------------
        if ((stagger_on) .and. g%periodic) then
            g%lu1i => g%memory(:, ig:)

            select case (g%mode_fdm1)
            case DEFAULT
                call FDM_C1INT6P_LHS(nx, g%jac, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3))
            end select
            call TRIDPFS(nx, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5))
            ig = ig + 5
        end if

! ###################################################################
! Check array sizes
! ###################################################################
        if (ig >= inb_grid) then
            call TLab_Write_ASCII(efile, __FILE__//'. Grid size incorrect.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'Leaving '//__FILE__)
#endif

        return
    end subroutine FDM_Initialize

end module FDM
