#include "dns_const.h"

!########################################################################
!#
!# Solving finite difference equations of the form
!#      (u')' - \lambda^2 u = f         ! regular case
!#      (u')' = f                       ! singular case
!#
!########################################################################
!# ARGUMENTS
!#
!# bcs     In    BCs data
!# u       Out   Solution
!# f       In    Forcing
!# tmp1    Out   First derivative
!#
!########################################################################
module OPR_FDE
    use TLAB_CONSTANTS, only: wp, wi
    implicit none
    private

    integer(wi) i
    real(wp) dummy, lambda
    real(wp), dimension(:), pointer :: a, b, c, d, e, g, h, ep, em

    integer, parameter :: i1 = 1, i2 = 2

    public :: FDE_BVP_SINGULAR_DD
    public :: FDE_BVP_SINGULAR_DN
    public :: FDE_BVP_SINGULAR_ND
    public :: FDE_BVP_SINGULAR_NN

    public :: FDE_BVP_REGULAR_DD
    public :: FDE_BVP_REGULAR_NN

contains
!########################################################################
!# Solving finite difference equations of the form
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!# See FDM_Int1_Initialize
!########################################################################
    subroutine OPR_Int1(g, f, result, ibc, bcs)
        use TLAB_CONSTANTS, only: BCS_MIN, BCS_MAX, BCS_BOTH
        use TLAB_TYPES, only: grid_dt
        use TLAB_ARRAYS, only: wrk2d
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: f(:, :)
        real(wp), intent(out) :: result(:, :)
        integer, intent(in) :: ibc
        real(wp), intent(in) :: bcs(:)

        integer(wi) :: idr, ndr, ic, len
        real(wp), pointer :: lhs(:, :) => null(), rhs(:, :) => null()

        ! ###################################################################
        len = size(f, 1)

        select case (ibc)
        case (BCS_MIN)
            result(:, 1) = bcs(:)
            result(:, g%size) = f(:, g%size)
            ! lhs =>
            ! rhs =>
        case (BCS_MAX)
            result(:, 1) = f(:, 1)
            result(:, g%size) = bcs(:)
            ! lhs =>
            ! rhs =>
        end select

        select case (g%nb_diag_1(1))
        case (3)
            ! call MatMul_3d(g%size, len, rhs(:, 1), rhs(:, 3), f, result, BCS_BOTH, rhs_b=g%rhsi_b(1:3, 0:3), rhs_t=g%rhsi_t(0:2, 1:4), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        case (5)
        end select

        select case (g%nb_diag_1(2))
        case (3)
            call TRIDSS(g%size - 2, len, lhs(2:, 1), lhs(2:, 2), lhs(2:, 3), result(:, 2:))
        case (5)
            call PENTADSS(g%size - 2, len, lhs(2:, 1), lhs(2:, 2), lhs(2:, 3), lhs(2:, 4), lhs(2:, 5), result(:, 2:))
        end select

        idr = g%nb_diag_1(2)/2 + 1
        ndr = g%nb_diag_1(2)

        if (any([BCS_MAX] == ibc)) then
            result(:, 1) = wrk2d(:, 1)
            do ic = 1, idr - 1
                result(:, 1) = result(:, 1) + lhs(1, idr + ic)*result(:, 1 + ic)
            end do
            result(:, 1) = result(:, 1) + lhs(1, 1)*result(:, 1 + ic)
        end if
        if (any([BCS_MIN] == ibc)) then
            result(:, g%size) = wrk2d(:, 2)
            do ic = 1, idr - 1
                result(:, g%size) = result(:, g%size) + lhs(g%size, idr - ic)*result(:, g%size - ic)
            end do
            result(:, g%size) = result(:, g%size) + lhs(g%size, ndr)*result(:, g%size - ic)
        end if

        return
    end subroutine OPR_INT1

!########################################################################
!Dirichlet/Neumann boundary conditions at imin/imax
!########################################################################
    subroutine FDE_BVP_SINGULAR_DN(imode_fdm, imax, jkmax, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, imax, jkmax
        real(wp), dimension(imax) :: dx
        real(wp), dimension(imax, 9), target :: wrk1d
        real(wp), dimension(jkmax, imax) :: u, f, tmp1
        real(wp), dimension(jkmax, 2) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 8)
        h => wrk1d(:, 9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_imax given
! -----------------------------------------------------------------------
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i2, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(imax - 1, a, b, c, d, e)

!   obtain v^0, array tmp1
            call PENTADSS(imax - 1, jkmax, a, b, c, d, e, tmp1)
            tmp1(:, imax) = 0.0_wp
            do i = 1, imax
                tmp1(:, i) = tmp1(:, i) + bcs(:, 2) ! add v_N to free array bcs(:,2)
            end do

!   obtain v^1, array wrk1d(:,6)
            call PENTADSS(imax - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i2, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(imax - 1, a, b, c, d, e, g, h)

!   obtain v^0, array tmp1
            call HEPTADSS(imax - 1, jkmax, a, b, c, d, e, g, h, tmp1)
            tmp1(:, imax) = 0.0_wp
            do i = 1, imax
                tmp1(:, i) = tmp1(:, i) + bcs(:, 2) ! add v_N to free array bcs(:,2)
            end do

!   obtain v^1, array wrk1d(:,6)
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        end if

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i1, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i1, dx, tmp1, u)
            call PENTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain u^0
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
            bcs(:, 2) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 2) = (bcs(:, 2) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, tmp1, u)
            call HEPTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain u^0
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
            bcs(:, 2) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 2) = (bcs(:, 2) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6M_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + d(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! Constraint
        dummy = 1.0_wp/(dummy - wrk1d(1, 6))
        bcs(:, 2) = (tmp1(:, 1) - bcs(:, 2))*dummy

! Result
        do i = 1, imax
            u(:, i) = u(:, i) + bcs(:, 2)*wrk1d(i, 7) + bcs(:, 1)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 2)*wrk1d(i, 6)
        end do

        return
    end subroutine FDE_BVP_SINGULAR_DN

!########################################################################
!Neumann/Dirichlet boundary conditions at imin/imax
!########################################################################
    subroutine FDE_BVP_SINGULAR_ND(imode_fdm, imax, jkmax, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, imax, jkmax
        real(wp), dimension(imax) :: dx
        real(wp), dimension(imax, 9), target :: wrk1d
        real(wp), dimension(jkmax, imax) :: u, f, tmp1
        real(wp), dimension(jkmax, 2) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 8)
        h => wrk1d(:, 9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
        f(:, imax) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i1, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i1, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(imax, 6) = dx(imax); wrk1d(imax - 1, 6) = e(imax)*dx(imax) ! for v^1
            call PENTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain v^0, array tmp1
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), tmp1(1, 2))
            tmp1(:, 1) = 0.0_wp
            do i = 1, imax
                tmp1(:, i) = tmp1(:, i) + bcs(:, 1) ! add v_1 to free array bcs(:,1)
            end do

!   obtain v^1, array wrk1d(:,6)
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 6))
            wrk1d(1, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(imax, 6) = dx(imax); wrk1d(imax - 1, 6) = g(imax)*dx(imax) ! for v^1
            call HEPTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain v^0, array tmp1
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), tmp1(1, 2))
            tmp1(:, 1) = 0.0_wp
            do i = 1, imax
                tmp1(:, i) = tmp1(:, i) + bcs(:, 1) ! add v_1 to free array bcs(:,1)
            end do

!   obtain v^1, array wrk1d(:,6)
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 6))
            wrk1d(1, 6) = 0.0_wp
        end if

! -----------------------------------------------------------------------
! solve for u in u' = v, u_N given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i2, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i2, dx, tmp1, u)
            call PENTADFS(imax - 1, a, b, c, d, e)

!   obtain u^0
            call PENTADSS(imax - 1, jkmax, a, b, c, d, e, u)
            bcs(:, 1) = u(:, imax); u(:, imax) = 0.0_wp
            bcs(:, 1) = (bcs(:, 1) + a(imax)*u(:, imax - 2) + b(imax)*u(:, imax - 1) + c(imax)*u(:, imax))/dx(imax) !u^(0)'_imax

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6_RHS(imax, i1, i2, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(imax - 1, i1, a, b, c, d, e, wrk1d(1, 7))
            dummy = wrk1d(imax, 7); wrk1d(imax, 7) = 0.0_wp
            dummy = (dummy + a(imax)*wrk1d(imax - 2, 7) + b(imax)*wrk1d(imax - 1, 7) + c(imax)*wrk1d(imax, 7))/dx(imax) ! u^(1)'_imax
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i2, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i2, dx, tmp1, u)
            call HEPTADFS(imax - 1, a, b, c, d, e, g, h)

!   obtain u^0
            call HEPTADSS(imax - 1, jkmax, a, b, c, d, e, g, h, u)
            bcs(:, 1) = u(:, imax); u(:, imax) = 0.0_wp
            bcs(:, 1) = (bcs(:, 1) + b(imax)*u(:, imax - 2) + c(imax)*u(:, imax - 1) + d(imax)*u(:, imax))/dx(imax) !u^(0)'_imax

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6M_RHS(imax, i1, i2, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 7))
            dummy = wrk1d(imax, 7); wrk1d(imax, 7) = 0.0_wp
            dummy = (dummy + b(imax)*wrk1d(imax - 2, 7) + c(imax)*wrk1d(imax - 1, 7) + d(imax)*wrk1d(imax, 7))/dx(imax) ! u^(1)'_imax
        end if

! Constraint
        dummy = 1.0_wp/(dummy - wrk1d(imax, 6))
        bcs(:, 1) = (tmp1(:, imax) - bcs(:, 1))*dummy

! Result
        do i = 1, imax
            u(:, i) = u(:, i) + bcs(:, 1)*wrk1d(i, 7) + bcs(:, 2)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 1)*wrk1d(i, 6)
        end do

        return
    end subroutine FDE_BVP_SINGULAR_ND

!########################################################################
!Dirichlet/Dirichlet boundary conditions at imin/imax
!########################################################################
    subroutine FDE_BVP_SINGULAR_DD(imode_fdm, imax, jkmax, x, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, imax, jkmax
        real(wp), dimension(imax) :: dx, x
        real(wp), dimension(imax, 9), target :: wrk1d
        real(wp), dimension(jkmax, imax) :: u, f, tmp1
        real(wp), dimension(jkmax, 3) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 8)
        h => wrk1d(:, 9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v = u' in (u')' = f , u'_imax given
! -----------------------------------------------------------------------
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i2, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(imax - 1, a, b, c, d, e)

!   obtain v^0, array tmp1
            call PENTADSS(imax - 1, jkmax, a, b, c, d, e, tmp1)
            tmp1(:, imax) = 0.0_wp

!   obtain v^1, array wrk1d(:,6)
            call PENTADSS(imax - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i2, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(imax - 1, a, b, c, d, e, g, h)

!     obtain v^0, array tmp1
            call HEPTADSS(imax - 1, jkmax, a, b, c, d, e, g, h, tmp1)
            tmp1(:, imax) = 0.0_wp

!     obtain v^1, array wrk1d(:,6)
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        end if
! -----------------------------------------------------------------------
! solve for u in u' v f, u_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i1, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i1, dx, tmp1, u)
            call PENTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain u^0
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, tmp1, u)
            call HEPTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain u^0
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6M_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! Constraint
        dummy = 1.0_wp/(dummy - wrk1d(1, 6))
        bcs(:, 3) = (tmp1(:, 1) - bcs(:, 3))*dummy

! BCs
        dummy = 1.0_wp/(x(imax) - x(1))
        bcs(:, 2) = (bcs(:, 2) - u(:, imax) - bcs(:, 3)*wrk1d(imax, 7) - bcs(:, 1))*dummy

! Result
        do i = 1, imax
            u(:, i) = u(:, i) + bcs(:, 3)*wrk1d(i, 7) + bcs(:, 2)*(x(i) - x(1)) + bcs(:, 1)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 3)*wrk1d(i, 6) + bcs(:, 2)
        end do

        return
    end subroutine FDE_BVP_SINGULAR_DD

!########################################################################
!Neumann/Neumann boundary conditions at imin/imax; must be compatible!
!########################################################################
    subroutine FDE_BVP_SINGULAR_NN(imode_fdm, imax, jkmax, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, imax, jkmax
        real(wp), dimension(imax) :: dx
        real(wp), dimension(imax, 7), target :: wrk1d
        real(wp), dimension(jkmax, imax) :: u, f, tmp1
        real(wp), dimension(jkmax, 2) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 6)
        h => wrk1d(:, 7)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(imax, i1, a, b, c, d, e)
            call INT_C1N6_RHS(imax, jkmax, i1, dx, f, tmp1)
            call PENTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:))
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), tmp1(1, 2))
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(imax, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, f, tmp1)
            call HEPTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), tmp1(1, 2))
        end if
        tmp1(:, 1) = 0.0_wp
        do i = 1, imax
            tmp1(:, i) = tmp1(:, i) + bcs(:, 1) ! this step assumes compatible problem
        end do
! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
!   same l.h.s. as before
            call INT_C1N6_RHS(imax, jkmax, i1, dx, tmp1, u)
!   same l.h.s. as before
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
!   same l.h.s. as before
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, tmp1, u)
!   same l.h.s. as before
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
        end if
        u(:, 1) = 0.0_wp ! this integration constant is free and set to zero

        return
    end subroutine FDE_BVP_SINGULAR_NN

!########################################################################
!Neumann/Neumann boundary conditions at imin/imax
!########################################################################
    subroutine FDE_BVP_REGULAR_NN(imode_fdm, imax, jkmax, cst, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, imax, jkmax
        real(wp) cst
        real(wp), dimension(imax) :: dx
        real(wp), dimension(imax, 12), target :: wrk1d
        real(wp), dimension(jkmax, imax) :: u, f, tmp1
        real(wp), dimension(jkmax, 3) :: bcs

! -----------------------------------------------------------------------
        real(wp) g_1, g_2, a_1, b_1, deti

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 11)
        h => wrk1d(:, 12)

        ep => wrk1d(:, 9)
        em => wrk1d(:, 10)

        lambda = SQRT(cst)

! #######################################################################
! -----------------------------------------------------------------------
! 1st step; solve for v^(0) and v^(1)
! -----------------------------------------------------------------------
        dummy = -lambda
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(imax, i2, dx, dummy, a, b, c, d, e, ep)
            call INT_C1N6_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(imax - 1, a, b, c, d, e)

!   obtain e^(+), array ep
            call PENTADSS(imax - 1, i1, a, b, c, d, e, ep)

!   obtain v^(0), array tmp1
            call PENTADSS(imax - 1, jkmax, a, b, c, d, e, tmp1)
            tmp1(:, imax) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call PENTADSS(imax - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(imax, i2, dx, dummy, a, b, c, d, e, g, h, ep)
            call INT_C1N6M_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(imax - 1, a, b, c, d, e, g, h)

!   obtain e^(+), array ep
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, ep)

!   obtain v^(0), array tmp1
            call HEPTADSS(imax - 1, jkmax, a, b, c, d, e, g, h, tmp1)
            tmp1(:, imax) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        end if
! -----------------------------------------------------------------------
! 2nd step; solve for u^(0) and u^(1) and u^(2:)
! -----------------------------------------------------------------------
        dummy = lambda
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(imax, i1, dx, dummy, a, b, c, d, e, em)
            call INT_C1N6_RHS(imax, jkmax, i1, dx, tmp1, u)
            call PENTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain e^(m), array em
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), em(2:))
            g_1 = (c(1)*em(1) + d(1)*em(2) + e(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6_RHS(imax, i1, i1, dx, ep, wrk1d(1, 8))
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + c(1)*wrk1d(1, 8) + d(1)*wrk1d(2, 8) + e(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))

!   BCs; intermediate step to save memory space
            dummy = 1.0_wp - wrk1d(imax, 8)*lambda
            bcs(:, 3) = tmp1(:, 1) - bcs(:, 1)
            bcs(:, 2) = u(:, imax)*lambda + bcs(:, 2)
            bcs(:, 1) = bcs(:, 2) + bcs(:, 3)*em(imax)! a_0 *det
            bcs(:, 2) = bcs(:, 2)*ep(1) + bcs(:, 3)*dummy   ! b_0 *det

            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(imax, i1, dx, dummy, a, b, c, d, e, g, h, em)
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, tmp1, u)
            call HEPTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain e^(m), array em
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), em(2:))
            g_1 = (d(1)*em(1) + e(1)*em(2) + g(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6M_RHS(imax, i1, i1, dx, ep, wrk1d(1, 8))
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + d(1)*wrk1d(1, 8) + e(1)*wrk1d(2, 8) + g(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))

!   BCs; intermediate step to save memory space
            dummy = 1.0_wp - wrk1d(imax, 8)*lambda
            bcs(:, 3) = tmp1(:, 1) - bcs(:, 1)
            bcs(:, 2) = u(:, imax)*lambda + bcs(:, 2)
            bcs(:, 1) = bcs(:, 2) + bcs(:, 3)*em(imax)! a_0 *det
            bcs(:, 2) = bcs(:, 2)*ep(1) + bcs(:, 3)*dummy   ! b_0 *det

            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6M_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + d(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! BCs; final step
        deti = 1.0_wp/(1.0_wp - ep(1)*em(imax) - wrk1d(imax, 8)*lambda) ! inverse of determinant
        bcs(:, 1) = bcs(:, 1)*deti ! a_0
        bcs(:, 2) = bcs(:, 2)*deti ! b_0
        a_1 = (lambda*wrk1d(imax, 7) + wrk1d(1, 6)*em(imax))*deti
        b_1 = (lambda*wrk1d(imax, 7)*ep(1) + wrk1d(1, 6)*(1.0_wp - wrk1d(imax, 8)*lambda))*deti

! Constraint
        dummy = 1.0_wp/(dummy + b_1*g_1 - wrk1d(1, 6) + a_1*g_2)
        bcs(:, 3) = (tmp1(:, 1) - bcs(:, 3) - bcs(:, 2)*g_1 - bcs(:, 1)*g_2)*dummy

        dummy = 1.0_wp/lambda
        bcs(:, 1) = bcs(:, 1) + bcs(:, 3)*a_1
        bcs(:, 2) = (bcs(:, 2) + bcs(:, 3)*b_1)*dummy

! Result
        do i = 1, imax
            u(:, i) = u(:, i) + bcs(:, 3)*wrk1d(i, 7) + bcs(:, 1)*wrk1d(i, 8) + bcs(:, 2)*em(i)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 3)*wrk1d(i, 6) + bcs(:, 1)*ep(i) - lambda*u(:, i)
        end do

        return
    end subroutine FDE_BVP_REGULAR_NN

!########################################################################
!Dirichlet/Dirichlet boundary conditions at imin/imax
!########################################################################
    subroutine FDE_BVP_REGULAR_DD(imode_fdm, imax, jkmax, cst, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, imax, jkmax
        real(wp) cst
        real(wp), dimension(imax) :: dx
        real(wp), dimension(imax, 12), target :: wrk1d
        real(wp), dimension(jkmax, imax) :: u, f, tmp1
        real(wp), dimension(jkmax, 3) :: bcs

! -----------------------------------------------------------------------
        real(wp) g_1, g_2, a_1, deti

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 11)
        h => wrk1d(:, 12)

        ep => wrk1d(:, 9)
        em => wrk1d(:, 10)

        lambda = SQRT(cst)

! #######################################################################
! -----------------------------------------------------------------------
! 1st step; solve for v^(0) and v^(1)
! -----------------------------------------------------------------------
        dummy = -lambda
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(imax, i2, dx, dummy, a, b, c, d, e, ep)
            call INT_C1N6_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(imax - 1, a, b, c, d, e)

!   obtain e^(+), array ep
            call PENTADSS(imax - 1, i1, a, b, c, d, e, ep)

!   obtain v^(0), array tmp1
            call PENTADSS(imax - 1, jkmax, a, b, c, d, e, tmp1)
            tmp1(:, imax) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call PENTADSS(imax - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(imax, i2, dx, dummy, a, b, c, d, e, g, h, ep)
            call INT_C1N6M_RHS(imax, jkmax, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(imax - 1, a, b, c, d, e, g, h)

!   obtain e^(+), array ep
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, ep)

!   obtain v^(0), array tmp1
            call HEPTADSS(imax - 1, jkmax, a, b, c, d, e, g, h, tmp1)
            tmp1(:, imax) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call HEPTADSS(imax - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(imax, 6) = 0.0_wp
        end if

! -----------------------------------------------------------------------
! 2nd step; solve for u^(0) and u^(1) and u^(2)
! -----------------------------------------------------------------------
        dummy = lambda
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(imax, i1, dx, dummy, a, b, c, d, e, em)
            call INT_C1N6_RHS(imax, jkmax, i1, dx, tmp1, u)
            call PENTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain e^(m), array em
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), em(2:))
            g_1 = (c(1)*em(1) + d(1)*em(2) + e(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6_RHS(imax, i1, i1, dx, ep, wrk1d(1, 8))
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + c(1)*wrk1d(1, 8) + d(1)*wrk1d(2, 8) + e(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call PENTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(imax, i1, dx, dummy, a, b, c, d, e, g, h, em)
            call INT_C1N6M_RHS(imax, jkmax, i1, dx, tmp1, u)
            call HEPTADFS(imax - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain e^(m), array em
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), em(2:))
            g_1 = (d(1)*em(1) + e(1)*em(2) + g(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6M_RHS(imax, i1, i1, dx, ep, wrk1d(1, 8))
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + d(1)*wrk1d(1, 8) + e(1)*wrk1d(2, 8) + g(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call HEPTADSS(imax - 1, jkmax, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6M_RHS(imax, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(imax - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + d(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! BCs; final step
        deti = 1.0_wp/wrk1d(imax, 8) ! inverse of determinant
        bcs(:, 2) = (bcs(:, 2) - u(:, imax) - bcs(:, 1)*em(imax))*deti ! a_0
! bcs(:,1) = bcs(:,1)                                       ! b_0/lambda
        a_1 = -wrk1d(imax, 7)*deti
! b_1 = 0.0_wp

! Constraint
        g_1 = g_1*lambda
        dummy = 1.0_wp/(dummy - wrk1d(1, 6) + a_1*g_2)
        bcs(:, 3) = (tmp1(:, 1) - bcs(:, 3) - bcs(:, 1)*g_1 - bcs(:, 2)*g_2)*dummy

        bcs(:, 2) = bcs(:, 2) + bcs(:, 3)*a_1 ! a = v_N
! bcs(:,1) = bcs(:,1)                ! b = u_1

! Result
        do i = 1, imax
            u(:, i) = u(:, i) + bcs(:, 3)*wrk1d(i, 7) + bcs(:, 2)*wrk1d(i, 8) + bcs(:, 1)*em(i)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 3)*wrk1d(i, 6) + bcs(:, 2)*ep(i) - lambda*u(:, i)
        end do

        return
    end subroutine FDE_BVP_REGULAR_DD

end module OPR_FDE
