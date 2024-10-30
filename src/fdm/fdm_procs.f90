#include "dns_error.h"

!########################################################################
! Building blocks to construct FDMs
! Based on Lagrange polynomial for non-uniform grids
! Calculation of RHS for different stencil lengths and bcs (periodic|biased)
!########################################################################
module FDM_PROCS
    use TLab_Constants
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public Pi                ! Product function defined over interval given by idx(:), Pi(x-x_j) for all j in idx
    public Pi_p              ! First-order derivative of Pi
    public Pi_pp_3           ! Second-order derivative when idx has only 3 points
    public Lag               ! Lagrange polynomials on idx(:) around i
    public Lag_p             ! First-order derivative of Lag
    public Lag_pp_3          ! Second-order derivative when idx has only 3 points

    public coef_e1n2_biased  ! coefficients for the biased, 2. order approximation to 1. order derivative
    public coef_e1n3_biased  ! coefficients for the biased, 3. order approximation to 1. order derivative

    public FDM_Bcs_Neumann      ! Initialize arrays to impose Neumann Bcs
    public FDM_Bcs_Reduce

! ###################################################################
! Compact parameters (1st derivative of 6th-order pentadiagonal); to be removed
! ###################################################################
    real(wp), public :: C1N6M_ALPHA, C1N6M_BETA
    real(wp), public :: C1N6M_ALPHA2, C1N6M_BETA2
    real(wp), public :: C1N6M_A, C1N6M_B, C1N6M_C
    real(wp), public :: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
    real(wp), public :: C1N6M_BD2, C1N6M_CD3

contains
    !########################################################################
    !########################################################################
    function Pi(x, j, idx) result(f)    ! Product function on interval idx(:) evaluated at x_j
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            f = f*(x(j) - x(idx(k)))
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Pi_p(x, j, idx) result(f)
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        real(wp) dummy
        integer(wi) k, m

        f = 0.0_wp
        do k = 1, size(idx)
            dummy = 1.0_wp
            do m = 1, size(idx)
                if (m /= k) then
                    dummy = dummy*(x(j) - x(idx(m)))
                end if
            end do
            f = f + dummy
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Pi_pp_3(x, j, idx) result(f)       ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        f = 2.0_wp*(x(j) - x(idx(1)) + x(j) - x(idx(2)) + x(j) - x(idx(3)))

        return
    end function

!########################################################################
!########################################################################
    function Lag(x, j, i, idx) result(f)        ! Lagrange polynomials on idx(:) around i evaluated at x_j
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f*(x(j) - x(idx(k)))/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

    ! -------------------------------------------------------------------
    function Lag_p(x, j, i, idx) result(f)        ! 1. derivative of Lagrange polynomials on idx(:) around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k, m
        real(wp) den, dummy

        den = 1.0_wp
        f = 0.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                dummy = 1.0_wp
                do m = 1, size(idx)
                    if (idx(m) /= i .and. m /= k) then
                        dummy = dummy*(x(j) - x(idx(m)))
                    end if
                end do
                f = f + dummy
                den = den*(x(i) - x(idx(k)))
            end if
        end do
        f = f/den

        return
    end function

    ! -------------------------------------------------------------------
    function Lag_pp_3(x, j, i, idx) result(f)    ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i, idx(:)
        real(wp) f

        integer(wi) k

        f = 2.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

!########################################################################
!########################################################################
    function coef_e1n3_biased(x, i, backwards) result(coef) ! first-order derivative, explicit, 3. order
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(4)

        integer(wi) stencil(4), k

        if (present(backwards)) then
            stencil = [i, i - 1, i - 2, i - 3]
        else
            stencil = [i, i + 1, i + 2, i + 3]
        end if

        do k = 1, size(stencil)
            coef(k) = Lag_p(x, i, stencil(k), stencil(:))
        end do
        ! if uniform, [ -11/6 3 -3/2 1/3 ]/h

        return
    end function

    ! -------------------------------------------------------------------
    function coef_e1n2_biased(x, i, backwards) result(coef) ! first-order derivative, explicit, 2. order
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(3)

        integer(wi) stencil(3), k

        if (present(backwards)) then
            stencil = [i, i - 1, i - 2]
        else
            stencil = [i, i + 1, i + 2]
        end if

        do k = 1, size(stencil)
            coef(k) = Lag_p(x, i, stencil(k), stencil(:))
        end do
        ! if uniform, [ -3/2 2 -1/2 ]/h

        return
    end function

! #######################################################################
! #######################################################################
    subroutine FDM_Bcs_Neumann(ibc, lhs, rhs, rhs_b, rhs_t)
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: lhs(:, :)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(inout) :: rhs_b(:, :), rhs_t(:, :)

        integer(wi) idl, ndl, idr, ndr, ir, ic, nx
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points

        ! For A_22, we need idl >= idr -1
        if (idl < idr - 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. LHS array is too small.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        ! For b_21, we need idr >= idl
        if (idr < idl) then
            call TLab_Write_ASCII(efile, __FILE__//'. RHS array is too small.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! -------------------------------------------------------------------
        if (any([BCS_ND, BCS_NN] == ibc)) then
            rhs_b(1:idr, 1:ndr) = rhs(1:idr, 1:ndr)

            dummy = 1.0_wp/rhs(1, idr)      ! normalize by r11

            ! reduced array B^R_{22}
            rhs_b(1, 1:ndr) = -rhs_b(1, 1:ndr)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = idr + 1, ndr        ! columns
                    rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) + rhs_b(1 + ir, idr - ir)*rhs_b(1, ic)
                end do
                ! longer stencil at the boundary
                ic = ndr + 1
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) + rhs_b(1 + ir, idr - ir)*rhs_b(1, 1)
            end do

            ! reduced array A^R_{22}
            lhs(1, 1:ndl) = lhs(1, 1:ndl)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) - rhs_b(1 + ir, idr - ir)*lhs(1, ic)
                end do
                ! vector a^R_{21} stored in rhs
                ic = idr
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir)*lhs(1, idl)
            end do

            ! finalize vector a^R_{21}
            do ir = 1, idl - 1
                ic = idr
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)
            end do

            ! store a_11/b_11 in rhs
            rhs_b(1, idr) = lhs(1, idl)

        end if

        if (any([BCS_DN, BCS_NN] == ibc)) then
            rhs_t(1:idr, 1:ndr) = rhs(nx - idr + 1:nx, 1:ndr)

            dummy = 1.0_wp/rhs(nx, idr)     ! normalize by rnn

            ! reduced array B^R_{11}
            rhs_t(idr, 1:ndr) = -rhs_t(idr, 1:ndr)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = 1, idr - 1          ! columns
                    rhs_t(idr - ir, ic + ir) = rhs(nx - ir, ic + ir) + rhs(nx - ir, idr + ir)*rhs_t(idr, ic)
                end do
                ! longer stencil at the boundary
                ic = 0
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir) + rhs(nx - ir, idr + ir)*rhs_t(idr, ndr)
            end do

            ! reduced array A^R_{11}
            lhs(nx, 1:ndl) = lhs(nx, 1:ndl)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = 1, idl - 1          ! columns
                    lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) - rhs(nx - ir, idr + ir)*lhs(nx, ic)
                end do
                ! vector a^R_{1n} stored in rhs
                ic = idr
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir)*lhs(nx, idl)
            end do

            ! finalize vector a^R_{1n}
            do ir = 1, idl - 1
                ic = idr
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir) - lhs(nx - ir, idl + ir)
            end do

            ! store a_nn/b_nn in rhs
            rhs_t(idr, idr) = lhs(nx, idl)

        end if

        return
    end subroutine FDM_Bcs_Neumann

! #######################################################################
! #######################################################################
    subroutine FDM_Bcs_Reduce(ibc, lhs, rhs, rhs_b, rhs_t)
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: lhs(:, :)
        real(wp), intent(in), optional :: rhs(:, :)
        real(wp), intent(inout), optional :: rhs_b(:, 0:), rhs_t(0:, :)

        integer(wi) idl, ndl, idr, ndr, ir, ic, nx, nx_t
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points
        nx_t = idr                      ! # grid points affected by bcs; for clarity

        ! -------------------------------------------------------------------
        if (any([BCS_MIN, BCS_BOTH] == ibc)) then
            dummy = 1.0_wp/lhs(1, idl)      ! normalize by l11

            ! reduced array A^R_{22}
            lhs(1, 1:ndl) = -lhs(1, 1:ndl)*dummy
            do ir = 1, idl - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lhs(1, ic)
                end do
                ic = ndl + 1                ! longer stencil at the boundary
                lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lhs(1, 1)
            end do

            ! reduced array B^R_{22}
            if (present(rhs_b)) then
                if (size(rhs_b, 1) < max(idl, idr + 1) .or. size(rhs_b, 2) < max(ndl, ndr)) then
                    call TLab_Write_ASCII(efile, __FILE__//'. rhs_b array is too small.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)
                end if

                rhs_b(1:max(idl, idr + 1), 1:ndr) = rhs(1:max(idl, idr + 1), 1:ndr)

                rhs_b(1, 1:ndr) = rhs_b(1, 1:ndr)*dummy
                do ir = 1, idl - 1              ! rows
                    do ic = idr, ndr            ! columns; ic = idr corresponds to vector b^R_{21}
                        rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)*rhs_b(1, ic)
                    end do
                    ic = ndr + 1                ! longer stencil at the boundary
                    rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)*rhs_b(1, 1)
                end do
            end if

        end if

        if (any([BCS_MAX, BCS_BOTH] == ibc)) then
            dummy = 1.0_wp/lhs(nx, idl)     ! normalize by lnn

            ! reduced array A^R_{11}
            lhs(nx, 1:ndl) = -lhs(nx, 1:ndl)*dummy
            do ir = 1, idl - 1              ! rows
                ic = 0                      ! longer stencil at the boundary
                lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lhs(nx, ndl)
                do ic = 1, idl - 1          ! columns
                    lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lhs(nx, ic)
                end do
            end do

            ! reduced array B^R_{11}
            if (present(rhs_t)) then
                if (size(rhs_t, 1) < max(idl, idr + 1) .or. size(rhs_t, 2) < max(ndl, ndr)) then
                    call TLab_Write_ASCII(efile, __FILE__//'. rhs_t array is too small.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)
                end if

                rhs_t(nx_t - max(idl, idr + 1) + 1:nx_t, 1:ndr) = rhs(nx - max(idl, idr + 1) + 1:nx, 1:ndr)

                rhs_t(nx_t, 1:ndr) = rhs_t(nx_t, 1:ndr)*dummy
                do ir = 1, idl - 1              ! rows
                    ic = 0                      ! columns; ic = 0 corresponds to longer stencil at the boundary
                    rhs_t(nx_t - ir, ic + ir) = rhs_t(nx_t - ir, ic + ir) - lhs(nx - ir, idl + ir)*rhs_t(nx_t, ndr)
                    do ic = 1, idr              ! ic = idr corresponds to vector b^R_{1n}
                        rhs_t(nx_t - ir, ic + ir) = rhs_t(nx_t - ir, ic + ir) - lhs(nx - ir, idl + ir)*rhs_t(nx_t, ic)
                    end do
                end do
            end if

        end if

        return
    end subroutine FDM_Bcs_Reduce

end module FDM_PROCS
