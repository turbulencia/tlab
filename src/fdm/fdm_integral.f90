#include "dns_error.h"

!########################################################################
! Initialize arrays to calculate boundary-value problems and integrals
! based on the compact schemes of the classes Au' = Bu and Au'' = Bu
! This probably should be OPR_ODE in operators, but need to disentangle it from global fdm_dt...
!########################################################################

module FDM_Integral
    use TLab_Constants, only: wp, wi, pi_wp, efile, BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use FDM_MatMul
    use FDM_PROCS
    implicit none
    private

    type, public :: fdm_integral_dt
        sequence
        integer mode_fdm1                           ! original finite-difference method
        real(wp), allocatable :: nodes(:)           ! nodes position; this and the previous information is redundant but useful in ODEs
        !                                             This info also allows to reproduce lhs and rhs, if needed
        !
        real(wp) :: lambda
        integer :: bc                               ! type of boundary condition, [ BCS_MIN, BCS_MAX ]
        real(wp) :: rhs_b(1:5, 0:7), rhs_t(0:4, 8)  ! # of diagonals is 7, # rows is 7/2+1
        real(wp), allocatable :: lhs(:, :)          ! Often overwritten to LU decomposition. Maybe add lu array and keep this as well.
        real(wp), allocatable :: rhs(:, :)
    end type fdm_integral_dt
    type(fdm_integral_dt), public :: fdm_Int0(2)    ! Global integral plan for lambda = 0

    public FDM_Int1_Initialize                      ! Prepare to solve u' +\lambda u = f
    public FDM_Int1_CreateSystem
    public FDM_Int1_Solve

    public FDM_Int2_Initialize                      ! Prepare to solve (u')' - \lamba^2 u = f
    public FDM_Int2_CreateSystem
    public FDM_Int2_Solve

contains
    !########################################################################
    !#
    !#     u'_i + \lamba u_i = f_i  N-1 eqns
    !#     u_1 or u_N given         1   eqn
    !#     Au' = Bu                 N   eqns
    !#
    !# starting from generic diagonal matrices A (lhs) and B (rhs).
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
    subroutine FDM_Int1_Initialize(x, lhs, rhs, lambda, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: lhs(:, :)               ! diagonals in lhs, or matrix A
        real(wp), intent(in) :: rhs(:, :)               ! diagonals in rhs, or matrix B
        real(wp), intent(in) :: lambda                  ! system constant
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, nd

        !########################################################################
        call FDM_Int1_CreateSystem(x, lhs, rhs, lambda, fdmi)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        nd = size(fdmi%lhs, 2)              ! # of diagonals

        select case (nd)
        case (3)
            call TRIDFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3))
        case (5)
            call PENTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5))
        case (7)
            call HEPTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7))
        end select

        return
    end subroutine FDM_Int1_Initialize

    !########################################################################
    !########################################################################
    subroutine FDM_Int1_CreateSystem(x, lhs, rhs, lambda, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: lhs(:, :)               ! diagonals in lhs, or matrix A
        real(wp), intent(in) :: rhs(:, :)               ! diagonals in rhs, or matrix B
        real(wp), intent(in) :: lambda                  ! system constant
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) i
        integer(wi) idl, ndl, idr, ndr, ir, nx
        real(wp) dummy, rhsr_b(5, 0:7), rhsr_t(0:4, 8)

        ! ###################################################################
        ndl = size(lhs, 2)
        idl = size(lhs, 2)/2 + 1        ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = size(rhs, 2)/2 + 1        ! center diagonal in rhs
        nx = size(lhs, 1)               ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%lambda = lambda

        if (allocated(fdmi%nodes)) deallocate (fdmi%nodes)
        allocate (fdmi%nodes(nx))
        fdmi%nodes(1:nx) = x(1:nx)

        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%lhs(nx, ndr))
        allocate (fdmi%rhs(nx, ndl))

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A), independent of lambda
        fdmi%rhs(:, :) = lhs(:, :)

        call FDM_Bcs_Reduce(fdmi%bc, fdmi%rhs, rhs, rhsr_b, rhsr_t)

        fdmi%rhs_b = 0.0_wp
        fdmi%rhs_t = 0.0_wp
        select case (fdmi%bc)
        case (BCS_MIN)
            fdmi%rhs_b(1:idl + 1, 1:ndl) = fdmi%rhs(1:idl + 1, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_b(1 + ir, idl - ir) = -rhsr_b(1 + ir, idr - ir)
            end do

        case (BCS_MAX)
            fdmi%rhs_t(0:idl, 1:ndl) = fdmi%rhs(nx - idl:nx, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_t(idl - ir, idl + ir) = -rhsr_t(idr - ir, idr + ir)
            end do

        end select

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C = B + h \lambda A), dependent on lambda
        fdmi%lhs(:, :) = rhs(:, :)

        fdmi%lhs(:, idr) = fdmi%lhs(:, idr) + lambda*lhs(:, idl)                ! center diagonal
        do i = 1, idl - 1                                                       ! off-diagonals
            fdmi%lhs(1 + i:nx, idr - i) = fdmi%lhs(1 + i:nx, idr - i) + lambda*lhs(1 + i:nx, idl - i)
            fdmi%lhs(1:nx - i, idr + i) = fdmi%lhs(1:nx - i, idr + i) + lambda*lhs(1:nx - i, idl + i)
        end do

        select case (fdmi%bc)
        case (BCS_MIN)
            fdmi%lhs(1:idr, 1:ndr) = rhsr_b(1:idr, 1:ndr)
            fdmi%lhs(1, idr + 1:idr + idl - 1) = fdmi%lhs(1, idr + 1:idr + idl - 1) - lambda*fdmi%rhs_b(1, idl + 1:ndl)
            ! fdmi%lhs(2:idr, 1:ndr) = rhsr_b(2:idr, 1:ndr)
            do ir = 1, idr - 1
                fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) + lambda*fdmi%rhs_b(1 + ir, 1:ndl)
            end do
        case (BCS_MAX)
            fdmi%lhs(nx - idr + 1:nx, 1:ndr) = rhsr_t(1:idr, 1:ndr)
            fdmi%lhs(nx, idr - idl + 1:idr - 1) = fdmi%lhs(nx, idr - idl + 1:idr - 1) - lambda*fdmi%rhs_t(idl, 1:idl - 1)
            ! fdmi%lhs(nx - idr + 1:nx - 1, 1:ndr) = rhsr_t(1:idr - 1, 1:ndr)
            do ir = 1, idr - 1
                fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) + lambda*fdmi%rhs_t(idl - ir, 1:ndl)
            end do
        end select

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 1, max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, idl)
            fdmi%rhs_b(ir, 0:ndl) = fdmi%rhs_b(ir, 0:ndl)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, idl)
            fdmi%rhs_t(idl - ir + 1, 1:ndl + 1) = fdmi%rhs_t(idl - ir + 1, 1:ndl + 1)*dummy

        end do

        ! do ir = 1, nx
        do ir = 2, nx - 1
            dummy = 1.0_wp/fdmi%rhs(ir, idl)

            fdmi%rhs(ir, 1:ndl) = fdmi%rhs(ir, 1:ndl)*dummy
            fdmi%lhs(ir, 1:ndr) = fdmi%lhs(ir, 1:ndr)*dummy

        end do

        ! -------------------------------------------------------------------
        ! reducing system in the opposite end to account for the case of extended stencils
        ! to move it up, you need to recalculate the expression for p_1 and p_n because they assume division by a_11
        select case (fdmi%bc)
        case (BCS_MIN)
            call FDM_Bcs_Reduce(BCS_MAX, fdmi%lhs, fdmi%rhs, rhs_t=fdmi%rhs_t)
        case (BCS_MAX)
            call FDM_Bcs_Reduce(BCS_MIN, fdmi%lhs, fdmi%rhs, rhs_b=fdmi%rhs_b)
        end select

        return
    end subroutine FDM_Int1_CreateSystem

    !########################################################################
    !########################################################################
    ! Allow to pass separate rhs because this part does not depend on lambda
    subroutine FDM_Int1_Solve(nlines, fdmi, rhsi, f, result, wrk2d, du_boundary)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: rhsi(:, :)
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2)
        real(wp), intent(out), optional :: du_boundary(nlines)

        ! -------------------------------------------------------------------
        integer(wi) :: nx
        integer(wi) :: idl, ndl, idr, ndr, ic

        ! ###################################################################
        nx = size(fdmi%lhs, 1)

        ndl = size(fdmi%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(rhsi, 2)
        idr = ndr/2 + 1

        select case (fdmi%bc)
        case (BCS_MIN)
            result(:, nx) = f(:, nx)
        case (BCS_MAX)
            result(:, 1) = f(:, 1)
        end select

        select case (size(rhsi, 2))
        case (3)
            call MatMul_3d(nx, nlines, rhsi(:, 1), rhsi(:, 3), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:3, 0:3), rhs_t=fdmi%rhs_t(0:2, 1:4), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        case (5)
            call MatMul_5d(nx, nlines, rhsi(:, 1), rhsi(:, 2), rhsi(:, 4), rhsi(:, 5), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:4, 0:5), rhs_t=fdmi%rhs_t(0:3, 1:6), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        end select

        select case (ndl)
        case (3)
            call TRIDSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), result(:, 2:))
        case (5)
            call PENTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), result(:, 2:))
        case (7)
            call HEPTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7), result(:, 2:))
        end select

        if (any([BCS_MAX] == fdmi%bc)) then
            result(:, 1) = wrk2d(:, 1) !*fdmi%lhs(1, idl)
            do ic = 1, idl - 1
                result(:, 1) = result(:, 1) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
            end do
            result(:, 1) = result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic)
            ! result(:, 1) = (result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic))*fdmi%lhs(1, idl)

            if (present(du_boundary)) then      ! calculate u'n
                du_boundary(:) = fdmi%lhs(nx, idl)*result(:, nx)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
                end do
                ic = idl                        ! longer stencil at the boundary
                du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, ndl)*result(:, nx - ic)

                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + rhsi(nx, idr - ic)*f(:, nx - ic)
                end do

            end if

        end if

        if (any([BCS_MIN] == fdmi%bc)) then
            result(:, nx) = wrk2d(:, 2) !*fdmi%lhs(nx, idl)
            do ic = 1, idl - 1
                result(:, nx) = result(:, nx) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
            end do
            result(:, nx) = result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic)
            ! result(:, nx) = (result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic))*fdmi%lhs(nx, idl)

            if (present(du_boundary)) then      ! calculate u'1
                du_boundary(:) = fdmi%lhs(1, idl)*result(:, 1)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
                end do
                ic = idl                        ! longer stencil at the boundary
                du_boundary(:) = du_boundary(:) + fdmi%lhs(1, 1)*result(:, 1 + ic)

                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + rhsi(1, idr + ic)*f(:, 1 + ic)
                end do

            end if

        end if

        return
    end subroutine FDM_Int1_Solve

    !########################################################################
    !#
    !#     u''_i - \lamba^2 u_i = f_i  N-2 eqns
    !#     u_1 and u_N given           2   eqns
    !#     Au'' = Bu                   N   eqns
    !#
    !# starting from the matrices A (lhs, tridiagonal) and B (rhs, pentadiagonal, with unitary central diagonal).
    !#
    !# The system of N-2 eqns:
    !#
    !#                    (B - \lambda^2 A)u = Af = g
    !#
    !# is established in this routine, giving diagonals a-e and g (see notes).
    !#
    !# The system is normalized such that the central diagonal in the new rhs is 1
    !#
    !########################################################################
    subroutine FDM_Int2_Initialize(x, lhs, rhs, lambda2, fdmi, f)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: lhs(:, :)               ! diagonals in lhs, or matrix A
        real(wp), intent(in) :: rhs(:, :)               ! diagonals in rhs, or matrix B
        real(wp), intent(in) :: lambda2                  ! system constant
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated
        real(wp), intent(out) :: f(:, :)                ! forcing terms for the hyperbolic sine

        ! -------------------------------------------------------------------
        integer(wi) nx, nd
        integer :: nlines = 1
        ! real(wp) bcs(2, 2)
        !########################################################################
        call FDM_Int2_CreateSystem(x, lhs, rhs, lambda2, fdmi, f)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        nd = size(fdmi%lhs, 2)              ! # of diagonals

        select case (nd)
        case (3)
            call TRIDFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3))
        case (5)
            ! We rely on this routines not changing a(2:3), b(2), e(ny-2:ny-1), d(ny-1)
            call PENTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5))
        case (7)
            call HEPTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7))
        end select

        ! Solve for particular solutions; to be rewritten in terms of _Int2_Solve
        select case (nd)
        case (3)
            call TRIDSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), f(2:, 1))
            call TRIDSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), f(2:, 2))
        case (5)
            call PENTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), f(2:, 1))
            call PENTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), f(2:, 2))
        case (7)
            call HEPTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7), f(2:, 1))
            call HEPTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7), f(2:, 2))
        end select

        ! ! Setting the RHS for the hyperbolic sine; the minus sign in included here to save ops
        ! f(:, :) = 0.0_wp; f(1, :) = 1.0_wp; f(nx, :) = 1.0_wp
        ! bcs(1, 1) = 1.0_wp; bcs(nx, 1) = 0.0_wp   ! s(-)
        ! bcs(1, 2) = 0.0_wp; bcs(nx, 2) = 1.0_wp   ! s(+)
        ! call FDM_Int2_Solve(2, fdmi, f, f, bcs, f)

        return
    end subroutine FDM_Int2_Initialize

    !########################################################################
    !########################################################################
    subroutine FDM_Int2_CreateSystem(x, lhs, rhs, lambda2, fdmi, si)
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: lhs(:, :)               ! diagonals in lhs, or matrix A
        real(wp), intent(in) :: rhs(:, :)               ! diagonals in rhs, or matrix B; center diagonal is not saved because it is 1
        real(wp), intent(in) :: lambda2                 ! system constant
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated
        real(wp), intent(out) :: si(:, :)               ! forcing terms for the hyperbolic sine

        ! -------------------------------------------------------------------
        integer(wi) i
        integer(wi) idl, ndl, idr, ndr, nx

        real(wp) dummy1, dummy2, pprime, coef(5), l2_inv, l2_min, l2_max

        ! ###################################################################
        ndl = size(lhs, 2)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = ndr/2 + 1             ! center diagonal in rhs; which is not saved because it is 1
        nx = size(lhs, 1)           ! # grid points

        ! check sizes
        if (ndr /= 5) then
            call TLab_Write_ASCII(efile, __FILE__//'. Only pentadiagonal case developed.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (ndl /= 3) then
            call TLab_Write_ASCII(efile, __FILE__//'. Only tridiagonal case developed.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%lambda = lambda2

        if (allocated(fdmi%nodes)) deallocate (fdmi%nodes)
        allocate (fdmi%nodes(nx))
        fdmi%nodes(1:nx) = x(1:nx)

        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%lhs(nx, ndr))
        allocate (fdmi%rhs(nx, ndl - 1))        ! Newcenter diagonal is also not saved because it is set to 1

        ! -------------------------------------------------------------------
        ! new lhs (array C22R); remember rhs center diagonal is not saved because it was 1
        fdmi%lhs(:, 1:idr - 1) = rhs(:, 1:idr - 1)
        fdmi%lhs(:, idr) = 1.0_wp
        fdmi%lhs(:, idr + 1:ndr) = rhs(:, idr:ndr - 1)

        fdmi%lhs(:, idr) = fdmi%lhs(:, idr) - lambda2*lhs(:, idl)               ! center diagonal
        do i = 1, idl - 1                                                       ! off-diagonals
            fdmi%lhs(:, idr - i) = fdmi%lhs(:, idr - i) - lambda2*lhs(:, idl - i)
            fdmi%lhs(:, idr + i) = fdmi%lhs(:, idr + i) - lambda2*lhs(:, idl + i)
        end do

        ! -------------------------------------------------------------------
        ! new tridiagonal rhs (array A22R); new central diagonal is 1 and I only need subdiagonal and superdiagonal
        fdmi%rhs(:, 1:idl - 1) = lhs(:, 1:idl - 1)
        fdmi%rhs(:, idl:ndl - 1) = lhs(:, idl + 1:ndl)

        ! Boundary corrections
        i = 2
        dummy1 = lhs(2, 1)/lhs(1, 2)
        fdmi%lhs(i, 3:5) = fdmi%lhs(i, 3:5) - dummy1*fdmi%lhs(i - 1, [4, 5, 1])
        l2_min = lhs(i, 2) - dummy1*lhs(1, 3)       ! central diagonal in reduced lhs
        fdmi%rhs(i, 1) = 0.0_wp                      ! See MatMul_3d

        i = nx - 1
        dummy2 = lhs(nx - 1, 3)/lhs(nx, 2)
        fdmi%lhs(i, 1:3) = fdmi%lhs(i, 1:3) - dummy2*fdmi%lhs(i + 1, [5, 1, 2])
        l2_max = lhs(i, 2) - dummy2*lhs(nx, 1)    ! central diagonal in reduced lhs
        fdmi%rhs(i, ndl - 1) = 0.0_wp                      ! See MatMul_3d

        ! Setting the RHS for the hyperbolic sine; the minus sign in included here to save ops
        si(:, :) = 0.0_wp

        si(1, 1) = 1.0_wp            ! b21R; this element is simply the solution at imin of s(-)
        si(2, 1) = -(rhs(2, 2) - dummy1)
        si(3, 1) = -rhs(3, 1)

        si(nx, 2) = 1.0_wp         ! b2nR; this element is simply the solution at nx of s(+)
        si(nx - 1, 2) = -(rhs(nx - 1, 3) - dummy2)
        si(nx - 2, 2) = -rhs(nx - 2, 4)

        ! -------------------------------------------------------------------
        ! Corrections to the BCS_DD to account for Neumann using third-order fdm for derivative at the boundary
        if (any([BCS_ND, BCS_NN] == fdmi%bc)) then
            ! Coefficients in FDM p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, 1)                ! second-order
            ! coef(1:4) = coef_e1n3_biased(x, 1)                ! third-order
            coef(1:5) = coef_c1n4_biased(x, 1)                ! fourth-order

            ! Data to calculate p_1 in terms of p_2, p_3, p_4 and p'_1 and p''_2
            fdmi%lhs(1, 2:5) = -coef([5, 2, 3, 4])/coef(1)              ! l_12 contains coefficient for p''_2
            fdmi%lhs(1, 3) = fdmi%lhs(1, 3) + lambda2*fdmi%lhs(1, 2)                ! d + lambda^2h^2 e in notes, e only 1 component
            pprime = 1.0_wp/coef(1)

            ! Derived coefficients; see notes
            fdmi%lhs(2, 3:5) = fdmi%lhs(2, 3:5) - si(2, 1)*fdmi%lhs(1, 3:5)          ! in reduced C matrix; the minus sign comes from def of f1
            fdmi%lhs(3, 2:4) = fdmi%lhs(3, 2:4) - si(3, 1)*fdmi%lhs(1, 3:5)

            l2_min = l2_min + fdmi%lhs(1, 2)*si(2, 1)                     ! in reduced A matrix; the plus sign comes from def of f1
            fdmi%rhs(3, 1) = fdmi%rhs(3, 1) + fdmi%lhs(1, 2)*si(3, 1)

            si(:, 1) = pprime*si(:, 1)                              ! for particular solutions

        end if
        if (any([BCS_DN, BCS_NN] == fdmi%bc)) then
            ! Coefficients in FDM p'_n = b_1 p_n + b_2 p_{n-1} + b_3 p_{n-2} +...
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, nx, backwards=.true.)
            ! coef(1:4) = coef_e1n3_biased(x, nx, backwards=.true.)
            coef(1:5) = coef_c1n4_biased(x, nx, backwards=.true.)

            ! Data to calculate p_n in terms of p_{n-1}, p_{n-2} and p'_n and p''_{n-1}
            fdmi%lhs(nx, [1, 2, 3, 5]) = -coef([4, 3, 2, 5])/coef(1)                    ! l_n5 contains coefficient for p''_{n-1}
            fdmi%lhs(nx, 3) = fdmi%lhs(nx, 3) + lambda2*fdmi%lhs(nx, 5)                         ! d + lambda^2h^2 e in notes, e only 1 component
            pprime = 1.0_wp/coef(1)

            ! Derived coefficients; see notes
            fdmi%lhs(nx - 1, 1:3) = fdmi%lhs(nx - 1, 1:3) - si(nx - 1, 2)*fdmi%lhs(nx, 1:3)      ! in reduced C matrix; the minus sign comes from def of f2
            fdmi%lhs(nx - 2, 2:4) = fdmi%lhs(nx - 2, 2:4) - si(nx - 2, 2)*fdmi%lhs(nx, 1:3)

            l2_max = l2_max + fdmi%lhs(nx, 5)*si(nx - 1, 2)                              ! in reduced A matrix; the plus sign comes from def of f2
            fdmi%rhs(nx - 2, ndl - 1) = fdmi%rhs(nx - 2, ndl - 1) + fdmi%lhs(nx, 5)*si(nx - 2, 2)

            si(:, 2) = pprime*si(:, 2)                                                    ! for particular solutions

        end if

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        i = 2
        l2_inv = 1.0_wp/l2_min

        fdmi%rhs(i, :) = fdmi%rhs(i, :)*l2_inv
        fdmi%lhs(i, :) = fdmi%lhs(i, :)*l2_inv
        si(i, :) = si(i, :)*l2_inv

        do i = 3, nx - 2
            l2_inv = 1.0_wp/lhs(i, 2)

            fdmi%rhs(i, :) = fdmi%rhs(i, :)*l2_inv
            fdmi%lhs(i, :) = fdmi%lhs(i, :)*l2_inv
            si(i, :) = si(i, :)*l2_inv

        end do

        i = nx - 1
        l2_inv = 1.0_wp/l2_max

        fdmi%rhs(i, :) = fdmi%rhs(i, :)*l2_inv
        fdmi%lhs(i, :) = fdmi%lhs(i, :)*l2_inv
        si(i, :) = si(i, :)*l2_inv

        return
    contains
        !########################################################################
        ! 1. derivatie of interpolation polynomial between equations (15) and (16)
        !    p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
        !
        ! Notation in Shukla and Zhong (2005), JCP, 204, 404–429 for the interpolation:
        !
        !       +                    I_n: set of points where the function and derivatives are given
        !   +---+---+---+---...
        !   +       +   +            I_m: set of points where only the function is given.
        !########################################################################
        function coef_c1n4_biased(x, i, backwards) result(coef)
            real(wp), intent(in) :: x(:)
            integer(wi), intent(in) :: i
            logical, intent(in), optional :: backwards
            real(wp) coef(5)

            real(wp) a2, b1, b2, b3, b4
            real(wp) dx1, dx3, dx4
            real(wp) D
            integer(wi) set_m(3), i1, i2, i3, i4

            i1 = i
            if (present(backwards)) then
                ! same as fowards, but changing the signs of the increments w.r.t. i
                ! To understand it, e.g., define a new variable k = -j, where k is the
                ! discrete variable moving around i
                i2 = i - 1
                i3 = i - 2
                i4 = i - 3
            else
                i2 = i + 1
                i3 = i + 2
                i4 = i + 3
            end if
            dx1 = x(i2) - x(i1)
            dx3 = x(i2) - x(i3)
            dx4 = x(i2) - x(i4)
            set_m = [i1, i3, i4]

            ! -------------------------------------------------------------------
            a2 = 0.5_wp*(Pi(x, i1, set_m) - dx1*Pi_p(x, i1, set_m))/Pi_p(x, i2, set_m)

            b2 = Pi_p(x, i1, set_m)*(2.0_wp*Pi_p(x, i2, set_m) + dx1*Pi_pp_3(x, i2, set_m)) &
                 - Pi(x, i1, set_m)*Pi_pp_3(x, i2, set_m)
            b2 = 0.5_wp*b2/Pi(x, i2, set_m)/Pi_p(x, i2, set_m)

            ! -------------------------------------------------------------------
            D = Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m)
            b1 = Lag(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + 2*dx1*Lag_p(x, i2, i1, set_m)) &
                 - dx1*Lag_p(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m))
            b1 = -b1/dx1/D

            D = Lag(x, i2, i3, set_m) + dx3*Lag_p(x, i2, i3, set_m)
            b3 = Lag(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + 2*dx1*Lag_p(x, i2, i3, set_m)) &
                 - dx1*Lag_p(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + dx1*Lag_p(x, i2, i3, set_m))
            b3 = -b3/dx3/D

            D = Lag(x, i2, i4, set_m) + dx4*Lag_p(x, i2, i4, set_m)
            b4 = Lag(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + 2*dx1*Lag_p(x, i2, i4, set_m)) &
                 - dx1*Lag_p(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + dx1*Lag_p(x, i2, i4, set_m))
            b4 = -b4/dx4/D

            coef = [b1, b2, b3, b4, a2]

            ! if uniform, we should have ( -29/6 54/6 -27/6 2/6 )/h and 3h
            ! print*, [b1, b2, b3, b4]*(x(2)-x(1))
            ! print*, a2/(x(2)-x(1))

            return
        end function

    end subroutine FDM_Int2_CreateSystem

    !########################################################################
    !########################################################################
    subroutine FDM_Int2_Solve(nlines, fdmi, si, f, bcs, result)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: si(:, :)                            ! particular solutions
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))

        ! -------------------------------------------------------------------
        integer(wi) :: nx, i
        integer(wi) :: ndl

        ! ###################################################################
        nx = size(fdmi%lhs, 1)

        ndl = size(fdmi%lhs, 2)

        ! Construct rhs
        ! note that current version does not retain the unitary diagonal,
        ! i.e. fdmi%rhs(:,2) is first upper diagonal in tridiagonal systems.
        result(:, 1) = 0.0_wp       ! This element is the solution at imin of p(0)
        result(:, nx) = 0.0_wp      ! This element is the solution at imax of p(0)
        call MatMul_3d(nx - 2, nlines, fdmi%rhs(2:, 1), fdmi%rhs(2:, 2), f(:, 2:), result(:, 2:))

        ! Solve pentadiagonal linear system
        select case (ndl)
        case (3)
            call TRIDSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), result(:, 2:))
        case (5)
            call PENTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), result(:, 2:))
        case (7)
            call HEPTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7), result(:, 2:))
        end select

        do i = 1, nx
            result(:, i) = result(:, i) + bcs(:, 1)*si(i, 1) + bcs(:, 2)*si(i, 2)
        end do

        !   Corrections to the BCS_DD to account for Neumann
        if (any([BCS_ND, BCS_NN] == fdmi%bc)) then
            result(:, 1) = result(:, 1) + fdmi%lhs(1, 3)*result(:, 2) &
                           + fdmi%lhs(1, 4)*result(:, 3) + fdmi%lhs(1, 5)*result(:, 4) &
                           + fdmi%lhs(1, 2)*f(:, 2)
        end if

        if (any([BCS_DN, BCS_NN] == fdmi%bc)) then
            result(:, nx) = result(:, nx) + fdmi%lhs(nx, 3)*result(:, nx - 1) &
                            + fdmi%lhs(nx, 2)*result(:, nx - 2) + fdmi%lhs(nx, 1)*result(:, nx - 3) &
                            + fdmi%lhs(nx, 5)*f(:, nx - 1)
        end if

        return
    end subroutine FDM_Int2_Solve

end module FDM_Integral
