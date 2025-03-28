#include "dns_error.h"

module FDM
    use TLab_Constants, only: wp, wi, roundoff_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
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
        integer(wi) i, ig, nx, ndl, ndr, inb_grid
        integer ib, bcs_cases(2)
        logical :: periodic_aux

        ! ###################################################################
        ! Consistency check
        ! ###################################################################
        if (g%periodic .and. g%mode_fdm1 == FDM_COM4_DIRECT) g%mode_fdm1 = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%mode_fdm1 == FDM_COM6_DIRECT) g%mode_fdm1 = FDM_COM6_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%der2%mode_fdm == FDM_COM4_DIRECT) g%der2%mode_fdm = FDM_COM4_JACOBIAN        ! they are the same for uniform grids.
        if (g%periodic .and. g%der2%mode_fdm == FDM_COM6_DIRECT) g%der2%mode_fdm = FDM_COM6_JACOBIAN_HYPER  ! they are the same for uniform grids.

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

        inb_grid = inb_grid + 2                 ! I need 2 aux array to calculate the Jacobian for 2. order derivative. To be removed.

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

        ig = ig + 1

        ig = ig + 2                     ! aux array to calculate the Jacobian for 2. order derivative. To be removed.

        ! ###################################################################
        ! first-order derivative
        ! ###################################################################
        g%lhs1 => g%memory(:, ig:)
        ig = ig + 5
        g%rhs1 => g%memory(:, ig:)
        ig = ig + 7
        if (g%periodic) then                                        ! modified wavenumber
            g%mwn1 => g%memory(:, ig)
            ig = ig + 1
        end if
        g%lu1 => g%memory(:, ig:)

        if (g%periodic) then
            ig = ig + g%nb_diag_1(1) + 2
        else                            ! biased,  different BCs
            ig = ig + 5*4
        end if

        ! -------------------------------------------------------------------
        ! uniform grid to calculate Jacobian (used for the stencils below and also as grid spacing in the code).
        g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
        g%jac(:, 1) = 1.0_wp

        call FDM_Der1_Initialize(g, periodic=.false., bcs_cases=[BCS_DD])

        ! Calculating derivative dxds into g%jac(:, 1)
        periodic_aux = g%periodic
        g%periodic = .false.
        call FDM_Der1_Solve(1, [0, 0], g, g%lu1, nodes, g%jac(:, 1), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...
        g%periodic = periodic_aux

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        g%nodes(:) = nodes(1:nx)

        call FDM_Der1_Initialize(g, periodic=g%periodic, bcs_cases=[BCS_DD, BCS_ND, BCS_DN, BCS_NN])

        if (g%periodic) g%mwn1(:) = g%mwn1(:)/g%jac(1, 1)           ! normalized by dx

        ! ###################################################################
        ! Check array sizes
        ! ###################################################################
        if (ig - 1 > inb_grid) then
            call TLab_Write_ASCII(efile, __FILE__//'. Grid size incorrect.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! ###################################################################
        ! second-order derivative
        ! ###################################################################
        ! -------------------------------------------------------------------
        ! uniform grid to calculate Jacobian (needed to set up the Jacobian formulations).
        g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
        g%jac(:, 2) = 1.0_wp
        g%jac(:, 3) = 0.0_wp

        call FDM_Der2_Initialize(g%nodes, g%jac(:, 2:), g%der2, periodic=.false., uniform=.true.)

        ! Calculating derivative d2xds2 into g%jac(:, 3)
        g%der2%periodic = .false.
        call FDM_Der2_Solve(1, g%der2, g%der2%lu, nodes, g%jac(:, 3), g%jac(:, 2), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        g%nodes(:) = nodes(1:nx)
        g%jac(:, 2) = g%jac(:, 1)

        call FDM_Der2_Initialize(g%nodes, g%jac(:, 2:), g%der2, g%periodic, g%uniform)

        if (g%der2%periodic) g%der2%mwn(:) = g%der2%mwn(:)/(g%jac(1, 1)**2)      ! normalized by dx

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
