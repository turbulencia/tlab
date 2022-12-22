#include "types.h"
#include "dns_const.h"

!########################################################################
!#
!# Calculate the boundary values of a field s.t. the normal derivative
!# is zero
!#
!# Routine format extracted from OPR_PARTIAL_Y
!#
!########################################################################
subroutine BOUNDARY_BCS_NEUMANN_Y(ibc, nx, ny, nz, g, u, bcs_hb, bcs_ht, wrk1d, tmp1, tmp2)
    use TLAB_TYPES, only: grid_dt
    use OPR_PARTIAL
    implicit none

    TINTEGER, intent(in) :: ibc     ! BCs at jmin/jmax: 1, for Neumann/-
    !                                                   2, for -      /Neumann
    !                                                   3, for Neumann/Neumann
    TINTEGER nx, ny, nz
    type(grid_dt), intent(IN) :: g
    TREAL, dimension(nx*nz, ny), target, intent(IN) :: u         ! they are transposed below
    TREAL, dimension(nx*nz, ny), target :: tmp1, tmp2 ! they are transposed below
    TREAL, dimension(g%size, 5), target :: wrk1d
    TREAL, dimension(nx*nz), target, intent(OUT) :: bcs_hb, bcs_ht

! -------------------------------------------------------------------
    TINTEGER nxz, nxy

    TREAL, dimension(:, :), pointer :: p_org, p_dst
    TREAL, dimension(:), pointer :: p_bcs_hb, p_bcs_ht
    TREAL, dimension(:), pointer :: a, b, c, d, e

! ###################################################################
    if (g%size == 1) then ! Set to zero in 2D case
        bcs_hb = C_0_R; bcs_ht = C_0_R

    else
! ###################################################################
        nxy = nx*ny
        nxz = nx*nz

        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)

! -------------------------------------------------------------------
! Make y  direction the last one
! -------------------------------------------------------------------
        if (nz == 1) then
            p_org => u
            p_dst => tmp1
            p_bcs_hb => bcs_hb
            p_bcs_ht => bcs_ht
        else
#ifdef USE_ESSL
            call DGETMO(u, nxy, nxy, nz, tmp1, nz)
#else
            call DNS_TRANSPOSE(u, nxy, nz, nxy, tmp1, nz)
#endif
            p_org => tmp1
            p_dst => tmp2
            p_bcs_hb => tmp1(:, 1)
            p_bcs_ht => tmp1(:, 2)
        end if

! ###################################################################
        select case (g%mode_fdm)

        case (FDM_COM4_JACOBIAN) !not yet implemented

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_BCS_LHS(ny, ibc, g%jac, a, b, c)
            call FDM_C1N6_BCS_RHS(ny, nxz, ibc, p_org, p_dst)

        case (FDM_COM6_JACPENTA)
            call FDM_C1N6M_BCS_LHS(ny, ibc, g%jac, a, b, c, d, e)
            call FDM_C1N6M_BCS_RHS(ny, nxz, ibc, p_org, p_dst)

        case (FDM_COM6_DIRECT) !not yet implemented
            call FDM_C1N6_BCS_LHS(ny, ibc, g%jac, a, b, c)
            call FDM_C1N6_BCS_RHS(ny, nxz, ibc, p_org, p_dst)

        case (FDM_COM8_JACOBIAN) !not yet implemented

        end select

! -------------------------------------------------------------------
        if (ibc == 1) then
            if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
                call TRIDFS(ny - 1, a(2:), b(2:), c(2:))
                call TRIDSS(ny - 1, nxz, a(2:), b(2:), c(2:), p_dst(:, 2))
                p_bcs_hb(:) = p_dst(:, 1) + c(1)*p_dst(:, 2)
            else
                call PENTADFS2(ny - 1, a(2:), b(2:), c(2:), d(2:), e(2:))
                call PENTADSS2(ny - 1, nxz, a(2:), b(2:), c(2:), d(2:), e(2:), p_dst(:, 2))
                p_bcs_hb(:) = p_dst(:, 1) + d(1)*p_dst(:, 2)
            end if

        else if (ibc == 2) then
            if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
                call TRIDFS(ny - 1, a, b, c)
                call TRIDSS(ny - 1, nxz, a, b, c, p_dst)
                p_bcs_ht(:) = p_dst(:, ny) + a(ny)*p_dst(:, ny - 1)
            else
                call PENTADFS2(ny - 1, a, b, c, d, e)
                call PENTADSS2(ny - 1, nxz, a, b, c, d, e, p_dst)
                p_bcs_ht(:) = p_dst(:, ny) + b(ny)*p_dst(:, ny - 1)
            end if

        else if (ibc == 3) then
            if (.not. (g%mode_fdm == FDM_COM6_JACPENTA)) then
                call TRIDFS(ny - 2, a(2:), b(2:), c(2:))
                call TRIDSS(ny - 2, nxz, a(2:), b(2:), c(2:), p_dst(:, 2))
                p_bcs_hb(:) = p_dst(:, 1) + c(1)*p_dst(:, 2)
                p_bcs_ht(:) = p_dst(:, ny) + a(ny)*p_dst(:, ny - 1)
            else
                call PENTADFS2(ny - 2, a(2:), b(2:), c(2:), d(2:), e(2:))
                call PENTADSS2(ny - 2, nxz, a(2:), b(2:), c(2:), d(2:), e(2:), p_dst(:, 2))
                p_bcs_hb(:) = p_dst(:, 1) + d(1)*p_dst(:, 2)
                p_bcs_ht(:) = p_dst(:, ny) + b(ny)*p_dst(:, ny - 1)
            end if
        end if

! ###################################################################
! -------------------------------------------------------------------
! Put bcs arrays in correct order
! -------------------------------------------------------------------
        if (nz > 1) then
#ifdef USE_ESSL
            if (ibc == 1 .or. ibc == 3) call DGETMO(p_bcs_hb, nz, nz, nx, bcs_hb, nx)
            if (ibc == 2 .or. ibc == 3) call DGETMO(p_bcs_ht, nz, nz, nx, bcs_ht, nx)
#else
            if (ibc == 1 .or. ibc == 3) call DNS_TRANSPOSE(p_bcs_hb, nz, nx, nz, bcs_hb, nx)
            if (ibc == 2 .or. ibc == 3) call DNS_TRANSPOSE(p_bcs_ht, nz, nx, nz, bcs_ht, nx)
#endif
        end if
        nullify (p_org, p_dst, p_bcs_hb, p_bcs_ht, a, b, c, d, e)

    end if

    return
end subroutine BOUNDARY_BCS_NEUMANN_Y
