#include "dns_error.h"

module FDM
    use TLab_Constants, only: wp, wi, roundoff_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    ! use TLab_Memory, only: TLab_Allocate_Real
    use FDM_PROCS, only: FDM_Bcs_Neumann
    use FDM_Derivative
    use FDM_Interpolate
    use FDM_Integral
    implicit none
    ! private                                   ! I need this to pass fdm_dt to parent modules
    !                                           to be fixed by decomposing fdm_dt into 3 derived types: der1, der2, interpol

    type(fdm_dt), public :: g(3)                ! Grid information along 3 directions

    public :: FDM_Initialize

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Initialize(nodes, g, fdmi, locScale)
        real(wp), intent(in) :: nodes(:)                            ! positions of the grid nodes
        type(fdm_dt), intent(inout) :: g                            ! fdm plan for derivatives
        type(fdm_integral_dt), intent(out), optional :: fdmi(2)     ! fdm plan for integrals
        real(wp), intent(in), optional :: locScale                  ! for consistency check

        ! -------------------------------------------------------------------
        integer(wi) i, ib, ip, ig, nx, ndl, ndr, inb_grid
        integer(wi) nmin, nmax, nsize, bcs_cases(4)

        integer, parameter :: i1 = 1
        logical :: periodic_aux

        ! ###################################################################
        ! Consistency check
        ! ###################################################################
        if (g%periodic .and. g%mode_fdm1 == FDM_COM4_DIRECT) g%mode_fdm1 = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm1 == FDM_COM6_DIRECT) g%mode_fdm1 = FDM_COM6_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm2 == FDM_COM4_DIRECT) g%mode_fdm2 = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm2 == FDM_COM6_DIRECT) g%mode_fdm2 = FDM_COM6_JACOBIAN_HYPER  ! they are the same for uniform grids.

        g%size = size(nodes)

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

        inb_grid = inb_grid + 1                 ! I need 1 aux array to calculate the Jacobian for 2. order derivative. To be removed.

        inb_grid = inb_grid &                   ! 1. order derivative
                   + 1 &                        ! Jacobian
                   + 5 &                        ! max # of diagonals in LHS
                   + 7                          ! max # of diagonals in RHS
        if (g%periodic) inb_grid = inb_grid + 1 ! modified wavenumber

        if (g%periodic) then
            inb_grid = inb_grid &
                       + 5 + 2                  ! LU decomposition 1. order
        else
            inb_grid = inb_grid &
                       + 5*4                    ! LU decomposition 1. order, 4 bcs
        end if

        inb_grid = inb_grid &                   ! 2. order derivative
                   + 1 &                        ! Jacobian
                   + 5 &                        ! max # of diagonals in LHS
                   + 7 + 5                      ! max # of diagonals in RHS
        if (g%periodic) inb_grid = inb_grid + 1 ! modified wavenumber

        if (g%periodic) then
            inb_grid = inb_grid &
                       + 5 + 2                  ! LU decomposition 2. order
        else
            inb_grid = inb_grid &
                       + 5                      ! LU decomposition 2. order, 1 bcs
        end if

        ! call TLab_Allocate_Real(__FILE__, g%memory, [g%size, inb_grid], g%name)
        allocate (g%memory(g%size, inb_grid))

        ! ###################################################################
        ! ###################################################################
        nx = g%size                     ! node number, for clarity below

        ig = 1                          ! Initialize counter to define pointers inside array x

        ! ###################################################################
        ! Space for node positions
        ! ###################################################################
        g%nodes => g%memory(:, ig)      ! Define pointer inside memory space

        ig = ig + 1                     ! Advance counter

        ! ###################################################################
        ! Space for Jacobians
        ! ###################################################################
        g%jac => g%memory(:, ig:)

        if (nx == 1) then
            g%jac(:, :) = 1.0_wp
            return
        end if

        ig = ig + 2

        ig = ig + 1         ! I need 1 aux array to calculate the Jacobian for 2. order derivative. To be removed.

        ! ###################################################################
        ! first-order derivative
        ! ###################################################################
        g%lhs1 => g%memory(:, ig:)
        ig = ig + 5
        g%rhs1 => g%memory(:, ig:)
        ig = ig + 7

        ! -------------------------------------------------------------------
        ! uniform grid to calculate Jacobian (used for the stencils below and also as grid spacing in the code).
        g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
        g%jac(:, 1) = 1.0_wp

        call FDM_Der1_CreateSystem(g, periodic=.false.)

        ! LU decomposition
        select case (g%nb_diag_1(1))
        case (3)
            call TRIDFS(nx, g%lhs1(:, 1), g%lhs1(:, 2), g%lhs1(:, 3))
        case (5)
            call PENTADFS2(nx, g%lhs1(:, 1), g%lhs1(:, 2), g%lhs1(:, 3), g%lhs1(:, 4), g%lhs1(:, 5))
        end select

        ! Calculating derivative dxds into g%jac(:, 1)
        periodic_aux = g%periodic
        g%periodic = .false.
        call FDM_Der1_Solve(1, [0, 0], g, g%lhs1, nodes, g%jac(:, 1), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...
        g%periodic = periodic_aux

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        g%nodes(:) = nodes(1:nx)

        if (g%periodic) then                                        ! modified wavenumber
            g%mwn1 => g%memory(:, ig)
            ig = ig + 1
        end if

        call FDM_Der1_CreateSystem(g, g%periodic)

        if (g%periodic) g%mwn1(:) = g%mwn1(:)/g%jac(1, 1)           ! normalized by dx

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

        ! ###################################################################
        ! second-order derivative
        ! ###################################################################
        g%lhs2 => g%memory(:, ig:)
        ig = ig + 5
        g%rhs2 => g%memory(:, ig:)
        ig = ig + 7 + 5

        ! -------------------------------------------------------------------
        ! uniform grid to calculate Jacobian (needed to set up the Jacobian formulations).
        g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
        g%jac(:, 2) = 1.0_wp
        g%jac(:, 3) = 0.0_wp

        call FDM_Der2_CreateSystem(g, periodic=.false.)

        ! LU decomposition
        select case (g%nb_diag_2(1))
        case (3)
            call TRIDFS(nx, g%lhs2(:, 1), g%lhs2(:, 2), g%lhs2(:, 3))
        case (5)
            call PENTADFS2(nx, g%lhs2(:, 1), g%lhs2(:, 2), g%lhs2(:, 3), g%lhs2(:, 4), g%lhs2(:, 5))
        end select

        ! Calculating derivative d2xds2 into g%jac(:, 3)
        periodic_aux = g%periodic
        g%periodic = .false.
        call FDM_Der2_Solve(1, g, g%lhs2, nodes, g%jac(:, 3), g%jac(:, 2), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...
        g%periodic = periodic_aux

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        g%nodes(:) = nodes(1:nx)
        g%jac(:, 2) = g%jac(:, 1)

        if (g%periodic) then                                        ! modified wavenumbers
            g%mwn2 => g%memory(:, ig)
            ig = ig + 1
        end if

        call FDM_Der2_CreateSystem(g, g%periodic)

        if (g%periodic) g%mwn2(:) = g%mwn2(:)/(g%jac(1, 1)**2)      ! normalized by dx

        ! -------------------------------------------------------------------
        ! LU decomposition
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

        ! ###################################################################
        ! Check array sizes
        ! ###################################################################
        if (ig >= inb_grid) then
            call TLab_Write_ASCII(efile, __FILE__//'. Grid size incorrect.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! ###################################################################
        ! interpolation for staggered cases
        ! ###################################################################
        if (stagger_on) then
            if (g%periodic) then

                call FDM_Interpol_Initialize(nodes(:), g%jac(:, 1), g%intl, g%mwn1(:))

                ! else
                !     call TLab_Write_ASCII(efile, 'Staggered grid only along periodic directions.')
                !     call TLab_Stop(DNS_ERROR_UNDEVELOP)

            end if

        end if

        ! ###################################################################
        ! first-order integrals (cases lambda = 0.0_wp)
        ! ###################################################################
        if (present(fdmi)) then
            if (g%periodic) then
                call TLab_Write_ASCII(wfile, __FILE__//'. Integral algorithms not available for periodic cases.')
            else
                ndl = g%nb_diag_1(1)
                ndr = g%nb_diag_1(2)

                bcs_cases(1:2) = [BCS_MIN, BCS_MAX]
                do ib = 1, 2
                    fdmi(ib)%mode_fdm = g%mode_fdm1
                    fdmi(ib)%bc = bcs_cases(ib)
                    call FDM_Int1_Initialize(g%nodes(:), g%lhs1(:, 1:ndl), g%rhs1(:, 1:ndr), 0.0_wp, fdmi(ib))

                end do
            end if
        end if

        return
    end subroutine FDM_Initialize

end module FDM
