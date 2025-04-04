!########################################################################
! Matrix multiplication of n-diagonal matrix with a vector
! The boundary conditions can extend over n/2 + 2 points
! The 1. upper-diagonal in interior points is equal to 1
! This allows use to handle systems A y = B x in which A amd B differ by up to 2 diagonals (see notes)
!########################################################################
module FDM_MatMul
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    implicit none
    private

    ! generic cases
    public MatMul_3d            ! Calculate f = B u, assuming B is tridiagonal
    public MatMul_3d_add        ! Calculate f = f + B u, assuming B is tridiagonal
    ! special cases where coefficients are constant in the interior points
    public MatMul_3d_antisym    ! Calculate f = B u, assuming B is tridiagonal, antisymmetric
    public MatMul_3d_sym        ! Calculate f = B u, assuming B is tridiagonal, symmetric

    ! generic cases
    public MatMul_5d            ! Calculate f = B u, assuming B is pentadiagonal
    public MatMul_5d_add        ! Calculate f = f + B u, assuming B is pentadiagonal
    ! special cases where coefficients are constant in the interior points
    public MatMul_5d_antisym    ! Calculate f = B u, assuming B is pentadiagonal, antisymmetric
    public MatMul_5d_sym        ! Calculate f = B u, assuming B is pentadiagonal, symmetric

    ! generic cases
    ! tbd when needed
    ! special cases where coefficients are constant in the interior points
    public MatMul_7d_antisym    ! Calculate f = B u, assuming B is heptadiagonal, antisymmetric
    public MatMul_7d_sym        ! Calculate f = B u, assuming B is heptadiagonal, symmetric

#define r0_b(j) rhs_b(j,0)
#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)
#define r4_b(j) rhs_b(j,4)
#define r5_b(j) rhs_b(j,5)
#define r6_b(j) rhs_b(j,6)
#define r7_b(j) rhs_b(j,7)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)
#define r5_t(j) rhs_t(j,5)
#define r6_t(j) rhs_t(j,6)
#define r7_t(j) rhs_t(j,7)

#define r1_i(j) rhs(j,1)
#define r2_i(j) rhs(j,2)
#define r3_i(j) rhs(j,3)
#define r4_i(j) rhs(j,4)
#define r5_i(j) rhs(j,5)
#define r6_i(j) rhs(j,6)
#define r7_i(j) rhs(j,7)

contains
    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is tri-diagonal and 1. upper-diagonal in interior points is equal to 1
    ! Special boundary conditions restricted to 3 points:
    ! r_11 r_12 r_13
    !      r_21 r_22 r_23
    !      r_30 r_31 r_32 r_33
    !                r_41  1.  r_43         <- interior points start here
    !                     ...  ...  ...
    subroutine MatMul_3d(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)  ! Special bcs at bottom and top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 3
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        ! -------------------------------------------------------------------
        ! Boundary; the first 3/2+1+1=3 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
            ! f(1) contains the boundary condition
            f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)
            f(:, 3) = f(:, 1)*r0_b(3) + u(:, 2)*r1_b(3) + u(:, 3)*r2_b(3) + u(:, 4)*r3_b(3)
        else
            f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r1_i(2) + u(:, 2)*r2_i(2) + u(:, 3)*r3_i(2)
            f(:, 3) = u(:, 2)*r1_i(3) + u(:, 3)*r2_i(3) + u(:, 4)*r3_i(3)
        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 4, nx - 3
            f(:, n) = u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n) + u(:, n + 1)
        end do

        ! -------------------------------------------------------------------
        ! Boundary; the last 3/2+1+1=3 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 2) = u(:, nx - 3)*r1_t(0) + u(:, nx - 2)*r2_t(0) + u(:, nx - 1)*r3_t(0) + f(:, nx)*r4_t(0)
            f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil
        else
            f(:, nx - 2) = u(:, nx - 3)*r1_i(nx - 2) + u(:, nx - 2)*r2_i(nx - 2) + u(:, nx - 1)*r3_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 2)*r1_i(nx - 1) + u(:, nx - 1)*r2_i(nx - 1) + u(:, nx)*r3_i(nx - 1)
            f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx) ! r3(nx) contains extended stencil
        end if

        return
    end subroutine MatMul_3d

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is tri-diagonal
    subroutine MatMul_3d_add(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! -------------------------------------------------------------------
        nx = size(rhs, 1)

        ! Boundary
        n = 1
        f(:, n) = f(:, n) + u(:, n)*r2_i(n) + u(:, n + 1)*r3_i(n) + u(:, n + 2)*r1_i(n)   ! r1(1) contains extended stencil

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = f(:, n) + u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n) + u(:, n + 1)*r3_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        f(:, n) = f(:, n) + u(:, n - 2)*r3_i(n) + u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n)   ! r3(nx) contains extended stencil

        return
    end subroutine MatMul_3d_add

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_antisym(rhs, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        integer ibc_loc

        ! -------------------------------------------------------------------
        nx = size(rhs, 1)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx)
            f(:, 2) = u(:, 3) - u(:, 1)

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains the boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
                f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)

            else
                f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1) ! r1(1) contains extended stencil
                f(:, 2) = u(:, 1)*r1_i(2) + u(:, 2)*r2_i(2) + u(:, 3)*r3_i(2)

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
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains the boundary condition
                f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
                if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil

            else
                f(:, nx - 1) = u(:, nx - 2)*r1_i(nx - 1) + u(:, nx - 1)*r2_i(nx - 1) + u(:, nx)*r3_i(nx - 1)
                f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx)  ! r3(nx) contains extended stencil

            end if

        end if

        return
    end subroutine MatMul_3d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_sym(rhs, u, f, periodic, ibc)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        integer ibc_loc

        ! -------------------------------------------------------------------
        nx = size(rhs, 1)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) + u(:, nx) + u(:, 1)*r2_i(1)

        else
            f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1) ! r1(1) contains extended stencil

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = u(:, n + 1) + u(:, n - 1) + u(:, n)*r2_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx) = u(:, 1) + u(:, nx - 1) + u(:, nx)*r2_i(nx)

        else
            f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx) ! r3(nx) contains extended stencil

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_3d_sym

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is penta-diagonal and 1. upper-diagonal in interior points is equal to 1
    subroutine MatMul_5d(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(1:4, 0:5), rhs_t(0:3, 1:6)  ! Special bcs at bottom and top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 5
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        ! -------------------------------------------------------------------
        ! Boundary; the first 5/2+1+1=4 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
            ! f(1) contains the boundary condition
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) contains extended stencil
            f(:, 2) = f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
            f(:, 3) = f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)
            f(:, 4) = f(:, 1)*r0_b(4) + u(:, 2)*r1_b(4) + u(:, 3)*r2_b(4) + u(:, 4)*r3_b(4) + u(:, 5)*r4_b(4) + u(:, 6)*r5_b(4)
        else
            f(:, 1) = u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) + u(:, 4)*r1_i(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)
            f(:, 3) = u(:, 1)*r1_i(3) + u(:, 2)*r2_i(3) + u(:, 3)*r3_i(3) + u(:, 4)*r4_i(3) + u(:, 5)*r5_i(3)
            f(:, 4) = u(:, 2)*r1_i(4) + u(:, 3)*r2_i(4) + u(:, 4)*r3_i(4) + u(:, 5)*r4_i(4) + u(:, 6)*r5_i(4)
        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 5, nx - 4
            f(:, n) = u(:, n - 2)*r1_i(n) + u(:, n - 1)*r2_i(n) + u(:, n)*r3_i(n) + u(:, n + 1) + u(:, n + 2)*r5_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary; the last 5/2+1+1=4 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 3) = u(:, nx - 5)*r1_t(0) + u(:, nx - 4)*r2_t(0) + u(:, nx - 3)*r3_t(0) + u(:, nx - 2)*r4_t(0) + u(:, nx - 1)*r5_t(0) + f(:, nx)*r6_t(0)
            f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
            f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) contains extended stencil
        else
        f(:, nx - 3) = u(:, nx - 5)*r1_i(nx - 3) + u(:, nx - 4)*r2_i(nx - 3) + u(:, nx - 3)*r3_i(nx - 3) + u(:, nx - 2)*r4_i(nx - 3) + u(:, nx - 1)*r5_i(nx - 3)
            f(:, nx - 2) = u(:, nx - 4)*r1_i(nx - 2) + u(:, nx - 3)*r2_i(nx - 2) + u(:, nx - 2)*r3_i(nx - 2) + u(:, nx - 1)*r4_i(nx - 2) + u(:, nx)*r5_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) + u(:, nx)*r4_i(nx - 1)
            f(:, nx) = u(:, nx - 3)*r5_i(nx) + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx) ! r5(nx) contains extended stencil
        end if

        return
    end subroutine MatMul_5d

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is pentadiagonal
    subroutine MatMul_5d_add(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! #######################################################################
        nx = size(rhs, 1)

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = f(:, 1) + u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) + u(:, 4)*r1_i(1)   ! r1(1) contains extended stencil
        f(:, 2) = f(:, 2) + u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = f(:, n) + u(:, n - 2)*r1_i(n) + u(:, n - 1)*r2_i(n) + u(:, n)*r3_i(n) + u(:, n + 1)*r4_i(n) + u(:, n + 2)*r5_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, nx - 1) = f(:, nx - 1) + u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) + u(:, nx)*r4_i(nx - 1)
        f(:, nx) = f(:, nx) + u(:, nx - 3)*r5_i(nx) + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx) ! r5(nx) contains extended stencil

        return
    end subroutine MatMul_5d_add

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. upper-diagonal for the interior points
    subroutine MatMul_5d_antisym(rhs, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r5_loc     ! 2. upper-diagonal
        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        r5_loc = r5_i(4)      ! The first 3 equations, last 3 equations, can be normalized differently

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) + r5_loc*(u(:, 3) - u(:, nx - 1))
            f(:, 2) = u(:, 3) - u(:, 1) + r5_loc*(u(:, 4) - u(:, nx))
            f(:, 3) = u(:, 4) - u(:, 2) + r5_loc*(u(:, 5) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) with extended stencil
                f(:, 2) = f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
                f(:, 3) = f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)

            else
                f(:, 1) = u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) + u(:, 4)*r1_i(1)   ! r1(1) with extended stencil
                f(:, 2) = u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)
                f(:, 3) = u(:, 1)*r1_i(3) + u(:, 2)*r2_i(3) + u(:, 3)*r3_i(3) + u(:, 4)*r4_i(3) + u(:, 5)*r5_i(3)

            end if

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r5_loc*(u(:, nx) - u(:, nx - 4))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r5_loc*(u(:, 1) - u(:, nx - 3))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r5_loc*(u(:, 2) - u(:, nx - 2))

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(n) contains boundary condition
                f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
                f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
                if (present(bcs_b)) bcs_t(:) = u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) with extended stencil

            else
            f(:, nx - 2) = u(:, nx - 4)*r1_i(nx - 2) + u(:, nx - 3)*r2_i(nx - 2) + u(:, nx - 2)*r3_i(nx - 2) + u(:, nx - 1)*r4_i(nx - 2) + u(:, nx)*r5_i(nx - 2)
                f(:, nx - 1) = u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) + u(:, nx)*r4_i(nx - 1)
                f(:, nx) = u(:, nx - 3)*r5_i(nx) + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx)! r5(nx) with extended stencil
            end if

        end if

        return
    end subroutine MatMul_5d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_5d_sym(rhs, u, f, periodic, ibc)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. upper-diagonal
        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        r5_loc = r5_i(3)      ! The first 2 equations, last 2 equations, are normalized differently
        r3_loc = r3_i(3)

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
            f(:, 1) = u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) &
                      + u(:, 4)*r1_i(1)   ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)

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
            f(:, nx - 1) = u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) &
                           + u(:, nx)*r4_i(nx - 1)

            f(:, nx) = u(:, nx - 3)*r5_i(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_5d_sym

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric heptadiagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. and 3. upper-diagonals for the interior points
    subroutine MatMul_7d_antisym(rhs, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r6_loc, r7_loc     ! 2. and 3. upper-diagonal
        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        r6_loc = r6_i(5)      ! The first 4 equations, last 4 equations, can be normalized differently
        r7_loc = r7_i(5)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) + r6_loc*(u(:, 3) - u(:, nx - 1)) + r7_loc*(u(:, 4) - u(:, nx - 2))
            f(:, 2) = u(:, 3) - u(:, 1) + r6_loc*(u(:, 4) - u(:, nx)) + r7_loc*(u(:, 5) - u(:, nx - 1))
            f(:, 3) = u(:, 4) - u(:, 2) + r6_loc*(u(:, 5) - u(:, 1)) + r7_loc*(u(:, 6) - u(:, nx))
            f(:, 4) = u(:, 5) - u(:, 3) + r6_loc*(u(:, 6) - u(:, 2)) + r7_loc*(u(:, 7) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r4_b(1) + u(:, 2)*r5_b(1) + u(:, 3)*r6_b(1) + u(:, 4)*r7_b(1) + u(:, 5)*r1_b(1)   ! r1(1) with extended stencil
                f(:, 2) = f(:, 1)*r3_b(2) + u(:, 2)*r4_b(2) + u(:, 3)*r5_b(2) + u(:, 4)*r6_b(2) + u(:, 5)*r7_b(2)
                f(:, 3) = f(:, 1)*r2_b(3) + u(:, 2)*r3_b(3) + u(:, 3)*r4_b(3) + u(:, 4)*r5_b(3) + u(:, 5)*r6_b(3) + u(:, 6)*r7_b(3)
                f(:, 4) = f(:, 1)*r1_b(4) + u(:, 2)*r2_b(4) + u(:, 3)*r3_b(4) + u(:, 4)*r4_b(4) + u(:, 5)*r5_b(4) + u(:, 6)*r6_b(4) + u(:, 7)*r7_b(4)

            else
                f(:, 1) = u(:, 1)*r4_i(1) + u(:, 2)*r5_i(1) + u(:, 3)*r6_i(1) + u(:, 4)*r7_i(1) + u(:, 5)*r1_i(1)   ! r1(1) with extended stencil
                f(:, 2) = u(:, 1)*r3_i(2) + u(:, 2)*r4_i(2) + u(:, 3)*r5_i(2) + u(:, 4)*r6_i(2) + u(:, 5)*r7_i(2)
                f(:, 3) = u(:, 1)*r2_i(3) + u(:, 2)*r3_i(3) + u(:, 3)*r4_i(3) + u(:, 4)*r5_i(3) + u(:, 5)*r6_i(3) + u(:, 6)*r7_i(3)
                f(:, 4) = u(:, 1)*r1_i(4) + u(:, 2)*r2_i(4) + u(:, 3)*r3_i(4) + u(:, 4)*r4_i(4) + u(:, 5)*r5_i(4) + u(:, 6)*r6_i(4) + u(:, 7)*r7_i(4)

            end if

        end if

        ! Interior points
        do n = 5, nx - 4
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r6_loc*(u(:, n + 2) - u(:, n - 2)) + r7_loc*(u(:, n + 3) - u(:, n - 3))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 3) = u(:, nx - 2) - u(:, nx - 4) + r6_loc*(u(:, nx - 1) - u(:, nx - 5)) + r7_loc*(u(:, nx) - u(:, nx - 6))
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r6_loc*(u(:, nx) - u(:, nx - 4)) + r7_loc*(u(:, 1) - u(:, nx - 5))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r6_loc*(u(:, 1) - u(:, nx - 3)) + r7_loc*(u(:, 2) - u(:, nx - 4))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r6_loc*(u(:, 2) - u(:, nx - 2)) + r7_loc*(u(:, 3) - u(:, nx - 3))

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains boundary condition
                f(:, nx - 3) = u(:, nx - 6)*r1_t(1) + u(:, nx - 5)*r2_t(1) + u(:, nx - 4)*r3_t(1) + u(:, nx - 3)*r4_t(1) + u(:, nx - 2)*r5_t(1) + u(:, nx - 1)*r6_t(1)+ f(:, nx)*r7_t(1)
              f(:, nx - 2) = u(:, nx - 5)*r1_t(2) + u(:, nx - 4)*r2_t(2) + u(:, nx - 3)*r3_t(2) + u(:, nx - 2)*r4_t(2) + u(:, nx - 1)*r5_t(2) + f(:, nx)*r6_t(2)
                f(:, nx - 1) = u(:, nx - 4)*r1_t(3) + u(:, nx - 3)*r2_t(3) + u(:, nx - 2)*r3_t(3) + u(:, nx - 1)*r4_t(3) + f(:, nx)*r5_t(3)
                if (present(bcs_b)) bcs_t(:) = u(:, nx - 4)*r7_t(4) + u(:, nx - 3)*r1_t(4) + u(:, nx - 2)*r2_t(4) + u(:, nx - 1)*r3_t(4) + f(:, nx)*r4_t(4) ! r7(nx) with extended stencil

            else
                f(:, nx - 3) = u(:, nx - 6)*r1_i(nx - 3) + u(:, nx - 5)*r2_i(nx - 3) + u(:, nx - 4)*r3_i(nx - 3) + u(:, nx - 3)*r4_i(nx - 3) + u(:, nx - 2)*r5_i(nx - 3) + u(:, nx - 1)*r6_i(nx - 3)+ u(:, nx)*r7_i(nx - 3)
                f(:, nx - 2) = u(:, nx - 5)*r1_i(nx - 2) + u(:, nx - 4)*r2_i(nx - 2) + u(:, nx - 3)*r3_i(nx - 2) + u(:, nx - 2)*r4_i(nx - 2) + u(:, nx - 1)*r5_i(nx - 2) + u(:, nx)*r6_i(nx - 2)
                f(:, nx - 1) = u(:, nx - 4)*r1_i(nx - 1) + u(:, nx - 3)*r2_i(nx - 1) + u(:, nx - 2)*r3_i(nx - 1) + u(:, nx - 1)*r4_i(nx - 1) + u(:, nx)*r5_i(nx - 1)
                f(:, nx) = u(:, nx - 4)*r7_i(nx) + u(:, nx - 3)*r1_i(nx) + u(:, nx - 2)*r2_i(nx) + u(:, nx - 1)*r3_i(nx) + u(:, nx)*r4_i(nx) ! r7(nx) with extended stencil
            end if

        end if

        return
    end subroutine MatMul_7d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_7d_sym(rhs, u, f, periodic, ibc)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal
        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        r7_loc = r7_i(4)
        r6_loc = r6_i(4)      ! The first 3 equations, last 3 equations, are normalized differently
        r4_loc = r4_i(4)

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
            f(:, 1) = u(:, 1)*r4_i(1) + u(:, 2)*r5_i(1) + u(:, 3)*r6_i(1) + u(:, 4)*r7_i(1) &
                      + u(:, 5)*r1_i(1)   ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r3_i(2) + u(:, 2)*r4_i(2) + u(:, 3)*r5_i(2) + u(:, 4)*r6_i(2) + u(:, 5)*r7_i(2)

            f(:, 3) = u(:, 1)*r2_i(3) + u(:, 2)*r3_i(3) + u(:, 3)*r4_i(3) + u(:, 4)*r5_i(3) + u(:, 5)*r6_i(3) + u(:, 6)*r7_i(3)

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
            f(:, nx - 2) = u(:, nx - 5)*r1_i(nx - 2) + u(:, nx - 4)*r2_i(nx - 2) + u(:, nx - 3)*r3_i(nx - 2) + u(:, nx - 2)*r4_i(nx - 2) + u(:, nx - 1)*r5_i(nx - 2) &
                           + u(:, nx)*r6_i(nx - 2)

            f(:, nx - 1) = u(:, nx - 4)*r1_i(nx - 1) + u(:, nx - 3)*r2_i(nx - 1) + u(:, nx - 2)*r3_i(nx - 1) + u(:, nx - 1)*r4_i(nx - 1) &
                           + u(:, nx)*r5_i(nx - 1)

            f(:, nx) = u(:, nx - 4)*r7_i(nx) & ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 3)*r1_i(nx) + u(:, nx - 2)*r2_i(nx) + u(:, nx - 1)*r3_i(nx) + u(:, nx)*r4_i(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_7d_sym

end module FDM_MatMul
