#include "dns_error.h"

!########################################################################
! Building blocks to construct FDMs
! Based on Lagrange polynomial for non-uniform grids
! Calculation of RHS for different stencil lengths and bcs (periodic|biased)
!########################################################################
module FDM_PROCS
    use TLAB_CONSTANTS
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
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

    ! generic cases
    public MatMul_3d            ! Calculate f = B u, assuming B is tridiagonal with center diagonal equal to 1
    public MatMul_3d_add        ! Calculate f = f + B u, assuming B is tridiagonal
    ! special cases, coefficients are constant in the interior points
    public MatMul_3d_antisym    ! Calculate f = B u, assuming B is tridiagonal, antisymmetric with 1. off-diagonal equal to 1
    public MatMul_3d_sym        ! Calculate f = B u, assuming B is tridiagonal, symmetric with 1. off-diagonal equal to 1
    !
    public MatMul_5d            ! Calculate f = B u, assuming B is pentadiagonal with center diagonal is 1
    public MatMul_5d_add        ! Calculate f = f + B u, assuming B is pentadiagonal
    public MatMul_5d_antisym    ! Calculate f = B u, assuming B is pentadiagonal, antisymmetric with 1. off-diagonal equal to 1
    public MatMul_5d_sym        ! Calculate f = B u, assuming B is pentadiagonal, symmetric with 1. off-diagonal equal to 1
    public MatMul_7d_sym        ! Calculate f = B u, assuming B is heptadiagonal, symmetric with 1. off-diagonal equal to 1
    public FDM_Bcs_Neumann
    public INT_C1NX_INITIALIZE

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
    function Pi(x, j, idx) result(f)    ! Product function defined over interval given by idx(:)
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
    function Lag(x, j, i, idx) result(f)        ! Lagrange polynomials on idx(:) around i
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
    ! Calculate f = B u, assuming B is tri-diagonal with center diagonal is 1
    subroutine MatMul_3d(nx, len, r1, r3, u, f)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r3(nx)  ! RHS diagonals (#=3-1 because center diagonal is 1)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1
        f(:, n) = u(:, n) + u(:, n + 1)*r3(n) &
                  + u(:, n + 2)*r1(n)   ! r11 contains 2. superdiagonal to allow for longer stencil at boundary

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = u(:, n - 1)*r1(n) + u(:, n) + u(:, n + 1)*r3(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        f(:, n) = u(:, n - 2)*r3(n) &   ! r3(n) contains 2. subdiagonal to allow for longer stencil at boundary
                  + u(:, n - 1)*r1(n) + u(:, n)

        return
    end subroutine MatMul_3d

    ! #######################################################################
    ! Calculate f = f + B u, assuming B is tri-diagonal
    subroutine MatMul_3d_add(nx, len, r1, r2, r3, u, f)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(inout) :: f(len, nx)    ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1
        f(:, n) = f(:, n) + u(:, n)*r2(n) + u(:, n + 1)*r3(n)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = f(:, n) + u(:, n - 1)*r1(n) + u(:, n)*r2(n) + u(:, n + 1)*r3(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        f(:, n) = f(:, n) + u(:, n - 1)*r1(n) + u(:, n)*r2(n)

        return
    end subroutine MatMul_3d_add

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric tri-diagonal with 1. off-diagonal equal to 1
    ! subroutine MatMul_3d_antisym(nx, len, r1, r2, r3, u, f, periodic, ibc)
    !     integer(wi), intent(in) :: nx, len                   ! m linear systems or size n
    !     real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)    ! RHS diagonals
    !     real(wp), intent(in) :: u(len, nx)                   ! function u
    !     real(wp), intent(out) :: f(len, nx)                  ! RHS, f = B u
    !     logical, intent(in) :: periodic
    !     integer, optional :: ibc

    !     ! -------------------------------------------------------------------
    !     integer(wi) n
    !     integer ibc_loc

    !     ! -------------------------------------------------------------------
    !     if (present(ibc)) then
    !         ibc_loc = ibc
    !     else
    !         ibc_loc = BCS_DD
    !     end if

    !     ! -------------------------------------------------------------------
    !     ! Boundary
    !     if (periodic) then
    !         f(:, 1) = u(:, 2) - u(:, nx)

    !     else
    !         f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) &
    !                   + u(:, 3)*r1(1)   ! r1(1) contains 2. superdiagonal to allow for longer stencil at boundary

    !         if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

    !     end if

    !     ! -------------------------------------------------------------------
    !     ! Interior points; accelerate
    !     do n = 2, nx - 1
    !         f(:, n) = u(:, n + 1) - u(:, n - 1)
    !     end do

    !     ! -------------------------------------------------------------------
    !     ! Boundary
    !     if (periodic) then
    !         f(:, nx) = u(:, 1) - u(:, nx - 1)

    !     else
    !         f(:, nx) = u(:, nx - 2)*r3(nx) & ! r3(nx) contains 2. subdiagonal to allow for longer stencil at boundary
    !                    + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx)

    !         if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

    !     end if

    !     return
    ! end subroutine MatMul_3d_antisym

    subroutine MatMul_3d_antisym(nx, len, r1, r2, r3, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len                   ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)    ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)                   ! function u
        real(wp), intent(inout) :: f(len, nx)                  ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc
        real(wp), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx)

            f(:, 2) = u(:, 3) - u(:, 1)

        else

            if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                if (present(bcs_b)) bcs_b(:) = u(:, 2)*r3_b(1) &
                                               + u(:, 3)*r1_b(1)   ! r1(1) contains 2. superdiagonal to allow for longer stencil at boundary

                f(:, 2) = u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2) + &
                          f(:, 1)*r1_b(2)                           ! f(1) contains the boundary condition
                ! bcs_b(:)*r1_b(2)

            else
                f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) &
                          + u(:, 3)*r1(1)   ! r1(1) contains 2. superdiagonal to allow for longer stencil at boundary

                f(:, 2) = u(:, 1)*r1(2) + u(:, 2)*r2(2) + u(:, 3)*r3(2)

            end if

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2)

            f(:, nx) = u(:, 1) - u(:, nx - 1)

        else
            if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + &
                               f(:, nx)*r3_t(1)

                if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) & ! r3(nx) contains 2. subdiagonal to allow for longer stencil at boundary
                                               + u(:, nx - 1)*r1_t(2)

            else
                f(:, nx - 1) = u(:, nx - 2)*r1(nx - 1) + u(:, nx - 1)*r2(nx - 1) + u(:, nx)*r3(nx - 1)

                f(:, nx) = u(:, nx - 2)*r3(nx) & ! r3(nx) contains 2. subdiagonal to allow for longer stencil at boundary
                           + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx)

            end if

        end if

#undef r1_b
#undef r2_b
#undef r3_b

#undef r1_t
#undef r2_t
#undef r3_t

        return
    end subroutine MatMul_3d_antisym

    ! #######################################################################
    ! Calculate f = B u, assuming B is symmetric tri-diagonal with 1. off-diagonal equal to 1
    subroutine MatMul_3d_sym(nx, len, r1, r2, r3, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, len                   ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)    ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)                   ! function u
        real(wp), intent(out) :: f(len, nx)                  ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) + u(:, nx) + u(:, 1)*r2(1)

        else
            f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) &
                      + u(:, 3)*r1(1)   ! r1(1) contains 2. superdiagonal to allow for longer stencil at boundary

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = u(:, n + 1) + u(:, n - 1) + u(:, n)*r2(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx) = u(:, 1) + u(:, nx - 1) + u(:, nx)*r2(nx)

        else
            f(:, nx) = u(:, nx - 2)*r3(nx) & ! r3(nx) contains 2. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_3d_sym

    ! #######################################################################
    ! Calculate f = B u, assuming B is penta-diagonal with center diagonal is 1
    subroutine MatMul_5d(nx, len, r1, r2, r3, r4, u, f)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx)    ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = u(:, 1) + u(:, 2)*r3(1) + u(:, 3)*r4(1) &
                  + u(:, 4)*r1(1)   ! r11 contains 3. superdiagonal to allow for longer stencil at boundary

        f(:, 2) = u(:, 1)*r2(2) + u(:, 2) + u(:, 3)*r3(2)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = u(:, n - 2)*r1(n) + u(:, n - 1)*r2(n) &
                      + u(:, n) &
                      + u(:, n + 1)*r3(n) + u(:, n + 2)*r4(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, nx - 1) = u(:, nx - 2)*r2(nx - 1) &
                       + u(:, nx - 1) &
                       + u(:, nx)*r3(nx - 1)

        f(:, nx) = u(:, nx - 3)*r4(nx) &   ! rn4 contains 3. subdiagonal to allow for longer stencil at boundary
                   + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)

        return
    end subroutine MatMul_5d

    ! #######################################################################
    ! Calculate f = f + B u, assuming B is pentadiagonal
    subroutine MatMul_5d_add(nx, len, r1, r2, r3, r4, r5, u, f)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(inout) :: f(len, nx)    ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = f(:, 1) + u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1)
        f(:, 2) = f(:, 2) + u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = f(:, n) + u(:, n - 2)*r1(n) + u(:, n - 1)*r2(n) + u(:, n)*r3(n) + u(:, n + 1)*r4(n) + u(:, n + 2)*r5(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, nx - 1) = f(:, nx - 1) + u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) + u(:, nx)*r4(nx - 1)
        f(:, nx) = f(:, nx) + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)

        return
    end subroutine MatMul_5d_add

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. off-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. off-diagonal for the interior points
    ! subroutine MatMul_5d_antisym(nx, len, r1, r2, r3, r4, r5, u, f, periodic, ibc)
    !     integer(wi), intent(in) :: nx, len       ! m linear systems or size n
    !     real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
    !     real(wp), intent(in) :: u(len, nx)       ! function u
    !     real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
    !     logical, intent(in) :: periodic
    !     integer, optional :: ibc

    !     ! -------------------------------------------------------------------
    !     integer(wi) n
    !     real(wp) r5_loc     ! 2. off-diagonal
    !     integer ibc_loc

    !     ! -------------------------------------------------------------------
    !     r5_loc = r5(3)      ! The first 2 equations, last 2 equations, are normalized differently

    !     if (present(ibc)) then
    !         ibc_loc = ibc
    !     else
    !         ibc_loc = BCS_DD
    !     end if

    !     ! Boundary
    !     if (periodic) then
    !         f(:, 1) = u(:, 2) - u(:, nx) &
    !                   + r5_loc*(u(:, 3) - u(:, nx - 1))

    !         f(:, 2) = u(:, 3) - u(:, 1) &
    !                   + r5_loc*(u(:, 4) - u(:, nx))

    !     else
    !         f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) &
    !                   + u(:, 4)*r1(1)   ! r1(1) contains 3. superdiagonal to allow for longer stencil at boundary

    !         f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

    !         if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

    !     end if

    !     ! Interior points
    !     do n = 3, nx - 2
    !         f(:, n) = u(:, n + 1) - u(:, n - 1) &
    !                   + r5_loc*(u(:, n + 2) - u(:, n - 2))
    !     end do

    !     ! Boundary
    !     if (periodic) then
    !         f(:, nx - 1) = u(:, nx) - u(:, nx - 2) &
    !                        + r5_loc*(u(:, 1) - u(:, nx - 3))

    !         f(:, nx) = u(:, 1) - u(:, nx - 1) &
    !                    + r5_loc*(u(:, 2) - u(:, nx - 2))

    !     else
    !         f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) &
    !                        + u(:, nx)*r4(nx - 1)
    !         f(:, nx) = u(:, nx - 3)*r5(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
    !                    + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)

    !         if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

    !     end if

    !     return
    ! end subroutine MatMul_5d_antisym

    subroutine MatMul_5d_antisym(nx, len, r1, r2, r3, r4, r5, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, len          ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)          ! function u
        real(wp), intent(inout) :: f(len, nx)       ! RHS, f = B u; f_1 and f_n can contain neumann bcs
        logical, intent(in) :: periodic
        integer, optional :: ibc
        real(wp), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), optional :: bcs_b(len), bcs_t(len)

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r5_loc     ! 2. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r5_loc = r5(4)      ! The first 2 equations, last 2 equations, are normalized differently

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)
#define r4_b(j) rhs_b(j,4)
#define r5_b(j) rhs_b(j,5)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)
#define r5_t(j) rhs_t(j,5)

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) &
                      + r5_loc*(u(:, 3) - u(:, nx - 1))

            f(:, 2) = u(:, 3) - u(:, 1) &
                      + r5_loc*(u(:, 4) - u(:, nx))

            f(:, 3) = u(:, 4) - u(:, 2) &
                      + r5_loc*(u(:, 5) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                if (present(bcs_b)) bcs_b(:) = u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1)     ! contribution to u_1

                f(:, 2) = u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2) + &
                          f(:, 1)*r2_b(2)       ! f(1) contains u'_1

                f(:, 3) = u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3) + &
                          f(:, 1)*r1_b(3)       ! f(1) contains u'_1

            else
                f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) &
                          + u(:, 4)*r1(1)   ! r1(1) contains 3. superdiagonal to allow for longer stencil at boundary

                f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

                f(:, 3) = u(:, 1)*r1(3) + u(:, 2)*r2(3) + u(:, 3)*r3(3) + u(:, 4)*r4(3) + u(:, 5)*r5(3)

            end if

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) &
                           + r5_loc*(u(:, nx) - u(:, nx - 4))

            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) &
                           + r5_loc*(u(:, 1) - u(:, nx - 3))

            f(:, nx) = u(:, 1) - u(:, nx - 1) &
                       + r5_loc*(u(:, 2) - u(:, nx - 2))

        else

            if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + &
                               f(:, nx)*r5_t(1)     ! f(n) contains u'_n

                f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + &
                               f(:, nx)*r4_t(2)     ! f(n) contains u'_n

                if (present(bcs_b)) bcs_t(:) = u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + u(:, nx - 3)*r5_t(3)    ! contribution to u_n

            else
                f(:, nx - 2) = u(:, nx - 4)*r1(nx - 2) + u(:, nx - 3)*r2(nx - 2) + u(:, nx - 2)*r3(nx - 2) &
                               + u(:, nx - 1)*r4(nx - 2) + u(:, nx)*r5(nx - 2)

                f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) &
                               + u(:, nx)*r4(nx - 1)

                f(:, nx) = u(:, nx - 3)*r5(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                           + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)
            end if

        end if

#undef r1_b
#undef r2_b
#undef r3_b
#undef r4_b
#undef r5_b

#undef r1_t
#undef r2_t
#undef r3_t
#undef r4_t
#undef r5_t

        return
    end subroutine MatMul_5d_antisym

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. off-diagonal equal to 1
    subroutine MatMul_5d_sym(nx, len, r1, r2, r3, r4, r5, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r5_loc = r5(3)      ! The first 2 equations, last 2 equations, are normalized differently
        r3_loc = r3(3)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = r3_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r5_loc*(u(:, 3) + u(:, nx - 1))

            f(:, 2) = r3_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r5_loc*(u(:, 4) + u(:, nx))

        else
            f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) &
                      + u(:, 4)*r1(1)   ! r1(1) contains 3. superdiagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 3, nx - 2
            f(:, n) = r3_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) + u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 1) = r3_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r5_loc*(u(:, 1) + u(:, nx - 3))

            f(:, nx) = r3_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r5_loc*(u(:, 2) + u(:, nx - 2))

        else
            f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) &
                           + u(:, nx)*r4(nx - 1)

            f(:, nx) = u(:, nx - 3)*r5(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_5d_sym

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric hepta-diagonal with 1. superdiagonal equal to 1
    subroutine MatMul_7d_sym(nx, len, r1, r2, r3, r4, r5, r6, r7, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, len       ! m linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx), r6(nx), r7(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(len, nx)       ! function u
        real(wp), intent(out) :: f(len, nx)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. off-diagonal
        real(wp) r7_loc     ! 3. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r7_loc = r7(4)
        r6_loc = r6(4)      ! The first 3 equations, last 3 equations, are normalized differently
        r4_loc = r4(4)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = r4_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r6_loc*(u(:, 3) + u(:, nx - 1)) &
                      + r7_loc*(u(:, 4) + u(:, nx - 2))

            f(:, 2) = r4_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r6_loc*(u(:, 4) + u(:, nx)) &
                      + r7_loc*(u(:, 5) + u(:, nx - 1))

            f(:, 3) = r4_loc*u(:, 3) + u(:, 4) + u(:, 2) &
                      + r6_loc*(u(:, 5) + u(:, 1)) &
                      + r7_loc*(u(:, 6) + u(:, nx))
        else
            f(:, 1) = u(:, 1)*r4(1) + u(:, 2)*r5(1) + u(:, 3)*r6(1) + u(:, 4)*r7(1) &
                      + u(:, 5)*r1(1)   ! r1(1) contains 4. superdiagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r3(2) + u(:, 2)*r4(2) + u(:, 3)*r5(2) + u(:, 4)*r6(2) + u(:, 5)*r7(2)

            f(:, 3) = u(:, 1)*r2(3) + u(:, 2)*r3(3) + u(:, 3)*r4(3) + u(:, 4)*r5(3) + u(:, 5)*r6(3) + u(:, 6)*r7(3)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = r4_loc*u(:, nx - 2) + u(:, nx - 1) + u(:, nx - 3) &
                           + r6_loc*(u(:, nx) + u(:, nx - 4)) &
                           + r7_loc*(u(:, 1) + u(:, nx - 5))

            f(:, nx - 1) = r4_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r6_loc*(u(:, 1) + u(:, nx - 3)) &
                           + r7_loc*(u(:, 2) + u(:, nx - 4))

            f(:, nx) = r4_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r6_loc*(u(:, 2) + u(:, nx - 2)) &
                       + r7_loc*(u(:, 3) + u(:, nx - 3))
        else
      f(:, nx - 2) = u(:, nx - 5)*r1(nx - 2) + u(:, nx - 4)*r2(nx - 2) + u(:, nx - 3)*r3(nx - 2) + u(:, nx - 2)*r4(nx - 2) + u(:, nx - 1)*r5(nx - 2) &
                           + u(:, nx)*r6(nx - 2)

            f(:, nx - 1) = u(:, nx - 4)*r1(nx - 1) + u(:, nx - 3)*r2(nx - 1) + u(:, nx - 2)*r3(nx - 1) + u(:, nx - 1)*r4(nx - 1) &
                           + u(:, nx)*r5(nx - 1)

            f(:, nx) = u(:, nx - 4)*r7(nx) & ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 3)*r1(nx) + u(:, nx - 2)*r2(nx) + u(:, nx - 1)*r3(nx) + u(:, nx)*r4(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_7d_sym

! ! #######################################################################
! ! #######################################################################
! ! TO BE IMPLEMENTED ACCORDING TO NOTES
!     subroutine FDM_Bcs(lhs, ibc)
!         real(wp), intent(inout) :: lhs(:, :)
!         integer, intent(in) :: ibc

!         integer(wi) id, nx

!         id = size(lhs, 2)/2 + 1     ! central diagonal
!         nx = size(lhs, 1)           ! # grid points

!         if (any([BCS_ND, BCS_NN] == ibc)) then
!             lhs(1, id:) = 0.0_wp
!             lhs(1, id) = 1.0_wp
!         end if

!         if (any([BCS_DN, BCS_NN] == ibc)) then
!             lhs(nx, :id) = 0.0_wp
!             lhs(nx, id) = 1.0_wp
!         end if

!         return
!     end subroutine FDM_Bcs

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

        if (idl < idr - 1) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. LHS is too small for Neumann BCs.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        ! -------------------------------------------------------------------
        if (any([BCS_ND, BCS_NN] == ibc)) then
            rhs_b(1:idr, 1:ndr) = rhs(1:idr, 1:ndr)

            dummy = 1.0_wp/rhs(1, idr)      ! normalize by r11

            rhs_b(1, 1:ndr) = -rhs(1, 1:ndr)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = idr + 1, ndr        ! columns
                    rhs_b(1 + ir, ic - ir) = rhs(1 + ir, ic - ir) + rhs(1 + ir, idr - ir)*rhs_b(1, ic)
                end do
                ! longer stencil at the boundary
                ic = ndr + 1
                rhs_b(1 + ir, ic - ir) = rhs_b(1 + ir, ic - ir) + rhs(1 + ir, idr - ir)*rhs_b(1, 1)
            end do

            lhs(1, :) = lhs(1, :)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    lhs(1 + ir, ic - ir) = lhs(1 + ir, ic - ir) - rhs(1 + ir, idr - ir)*lhs(1, ic)
                end do
                ! term for nonzero derivative
                rhs_b(1 + ir, idr - ir) = rhs_b(1 + ir, idr - ir)*lhs(1, idl)
            end do

            ! finalize term for nonzero derivative
            do ir = 1, idl - 1
                rhs_b(1 + ir, idr - ir) = rhs_b(1 + ir, idr - ir) - lhs(1 + ir, idl - ir)
            end do

        end if

        if (any([BCS_DN, BCS_NN] == ibc)) then
            nx = size(lhs, 1)               ! # grid points
            rhs_t(1:idr, 1:ndr) = rhs(nx - idr + 1:nx, 1:ndr)

            dummy = 1.0_wp/rhs(nx, idr)     ! normalize by rnn

            rhs_t(idr, :) = -rhs(nx, :)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = 1, idr - 1          ! columns
                    rhs_t(idr - ir, ic + ir) = rhs(nx - ir, ic + ir) + rhs(nx - ir, idr + ir)*rhs_t(idr, ic)
                end do
                ! longer stencil at the boundary
                ic = 0
                rhs_t(idr - ir, ic + ir) = rhs_t(idr - ir, ic + ir) + rhs(nx - ir, idr + ir)*rhs_t(idr, ndr)
            end do

            lhs(nx, :) = lhs(nx, :)*dummy
            do ir = 1, idr - 1              ! rows
                do ic = 1, idl - 1          ! columns
                    lhs(nx - ir, ic + ir) = lhs(nx - ir, ic + ir) - rhs(nx - ir, idr + ir)*lhs(nx, ic)
                end do
                ! term for nonzero derivative
                rhs_t(idr - ir, idr + ir) = rhs_t(idr - ir, idr + ir)*lhs(nx, idl)
            end do

            ! finalize term for nonzero derivative
            do ir = 1, idl - 1
                rhs_t(idr - ir, idr + ir) = rhs_t(idr - ir, idr + ir) - lhs(nx - ir, idl + ir)
            end do

        end if

        return
    end subroutine FDM_Bcs_Neumann

    !########################################################################
!#
!# Initialize the solver for the BVP
!#
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# The system of N-1 eqns:
!#
!#                    (B + \lambda A)u = Af
!#
!# is established in this routine (see notes).
!#
!# The system is normalized such that the central diagonal in the new rhs is 1
!#
!########################################################################
    subroutine INT_C1NX_INITIALIZE(ibc, lhs, rhs, lambda, lu, rhs_int)
        use TLAB_CONSTANTS
        implicit none

        integer, intent(in) :: ibc              ! Boundary condition, BCS_DN, BCS_ND
        real(wp), intent(in) :: lhs(:, :)       ! diagonals in lhs, or matrix A
        real(wp), intent(in) :: rhs(:, :)       ! diagonals in rhs, or matrix B
        real(wp), intent(in) :: lambda          ! system constant
        real(wp), intent(out) :: lu(:, :)       ! diagonals in new lhs
        real(wp), intent(out) :: rhs_int(:, :)  ! diagonals in new rhs

! -------------------------------------------------------------------
        integer(wi) i
        integer(wi) idl, ndl, idr, ndr, ir, ic, nx
        real(wp) dummy

        ! -------------------------------------------------------------------
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points

! ###################################################################
        ! check sizes
        if (size(lu, 2) < ndr) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong array lu size.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (size(rhs_int, 2) < ndl) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong array rhs_int size.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (idr < idl) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. New LHS is too small for integral operator.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        ! new lhs diagonals (array C22R)
        lu(:, 1:ndr) = rhs(:, 1:ndr)

        lu(:, idr) = lu(:, idr) + lambda*lhs(:, idl)    ! center diagonal
        do i = 1, idl - 1                               ! off-diagonals
            lu(1 + i:nx, idr - i) = lu(1 + i:nx, idr - i) + lambda*lhs(1 + i:nx, idl - i)      ! skip the top-left corner
            lu(1:nx - i, idr + i) = lu(1:nx - i, idr + i) + lambda*lhs(1:nx - i, idl + i)      ! skip the bottom-right corner
        end do

        ! new rhs diagonals (array A22R)
        rhs_int(:, 1:ndl) = lhs(:, 1:ndl)
        rhs_int(1, 1) = 0.0_wp              ! longer stencils at the boundaries
        rhs_int(nx, ndl) = 0.0_wp

        ! Boundary corrections
        select case (ibc)
        case (BCS_DN)
            dummy = 1.0_wp/lhs(1, idl)      ! normalize by l11

            ! term for nonzero bc
            do ir = 1, idr - 1
                lu(1 + ir, idr - ir) = -rhs(1 + ir, idr - ir)
            end do

            ! reduced array C22R
            lu(1, 1:ndr) = -lu(1, 1:ndr)*dummy
            do ir = 1, idl - 1              ! rows
                do ic = idr + 1, ndr        ! columns
                    lu(1 + ir, ic - ir) = lu(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lu(1, ic)
                end do
                ! longer stencil at the boundary
                ic = ndr + 1
                lu(1 + ir, ic - ir) = lu(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*lu(1, 1)
                ! term for nonzero bc
                ic = idr
                lu(1 + ir, ic - ir) = lu(1 + ir, ic - ir) + lhs(1 + ir, idl - ir)*rhs(1, idr)*dummy
            end do

            ! reduced array A22R
            rhs_int(1, 1:ndl) = rhs_int(1, 1:ndl)*dummy
            do ir = 1, idl - 1              ! rows
                do ic = idl + 1, ndl        ! columns
                    rhs_int(1 + ir, ic - ir) = rhs_int(1 + ir, ic - ir) - lhs(1 + ir, idl - ir)*rhs_int(1, ic)
                end do
            end do
            ! longer stencil at the boundary
            rhs_int(2, 1) = 0.0_wp

        case (BCS_ND)
            dummy = 1.0_wp/lhs(nx, idl)     ! normalize by lnn

            ! term for nonzero bc
            do ir = 1, idr - 1
                lu(nx - ir, idr + ir) = -rhs(nx - ir, idr + ir)
            end do

            ! reduced array C22R
            lu(nx, 1:ndr) = -lu(nx, 1:ndr)*dummy
            do ir = 1, idl - 1              ! rows
                do ic = 1, idr - 1          ! columns
                    lu(nx - ir, ic + ir) = lu(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lu(nx, ic)
                end do
                ! longer stencil at the boundary
                ic = 0
                lu(nx - ir, ic + ir) = lu(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*lu(nx, ndr)
                ! term for nonzero bc
                ic = idr
                lu(nx - ir, ic + ir) = lu(nx - ir, ic + ir) + lhs(nx - ir, idl + ir)*rhs(nx, idr)*dummy
            end do

            ! reduced array A22R
            rhs_int(nx, 1:ndl) = rhs_int(nx, 1:ndl)*dummy
            do ir = 1, idl - 1              ! rows
                do ic = 1, idl - 1          ! columns
                    rhs_int(nx - ir, ic + ir) = rhs_int(nx - ir, ic + ir) - lhs(nx - ir, idl + ir)*rhs_int(nx, ic)
                end do
            end do
            ! longer stencil at the boundary
            rhs_int(nx - 1, ndl) = 0.0_wp

        end select

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 1, nx
            dummy = 1.0_wp/rhs_int(ir, idl)

            rhs_int(ir, 1:ndl) = rhs_int(ir, 1:ndl)*dummy
            lu(ir, 1:ndr) = lu(ir, 1:ndr)*dummy

        end do

        return
    end subroutine INT_C1NX_INITIALIZE

end module FDM_PROCS
