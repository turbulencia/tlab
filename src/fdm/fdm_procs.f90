!########################################################################
! Building blocks to construct various FDM methods
! Based on Lagrange polynomial for non-uniform grids, works of course as well when uniform grids
! Calculation of RHS for different stencil lengths and bcs (periodic|biased)
!########################################################################
module FDM_PROCS
    use TLAB_CONSTANTS
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

    public MatMul_5d    ! Calculate f = B u, assuming B is pentadiagonal with center diagonal is 1
    public MatMul_3d      ! Calculate f = B u, assuming B is tridiagonal with center diagonal is 1
    public MatMul_5d_antisym
    public MatMul_5d_sym
    public MatMul_7d_sym
    public FDM_Bcs

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
    ! Calculate f = B u, assuming B is penta-diagonal with center diagonal is 1
    subroutine MatMul_3d(nmax, mmax, r1, r2, u, f)
        integer(wi), intent(in) :: nmax, mmax       ! m linear systems or size n
        real(wp), intent(in) :: r1(nmax), r2(nmax)  ! RHS diagonals (#=3-1 because center diagonal is 1)
        real(wp), intent(in) :: u(mmax, nmax)       ! function u
        real(wp), intent(out) :: f(mmax, nmax)      ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1 ! rhs(1,1) contains 1. superdiagonal to allow for longer stencil at boundary
        f(:, n) = &
            +u(:, n) &
            + u(:, n + 1)*r2(n) + u(:, n + 2)*r1(n)

        ! Interior points
        do n = 2, nmax - 1
            f(:, n) = u(:, n - 1)*r1(n) &
                      + u(:, n) &
                      + u(:, n + 1)*r2(n)
        end do

        ! Boundary
        n = nmax ! rhs(n,2) contains 1. subdiagonal to allow for longer stencil at boundary
        f(:, n) = u(:, n - 2)*r2(n) + u(:, n - 1)*r1(n) &
                  + u(:, n)

        return
    end subroutine MatMul_3d

    ! #######################################################################
    ! Calculate f = B u, assuming B is penta-diagonal with center diagonal is 1
    subroutine MatMul_5d(nmax, mmax, rhs, u, f)
        integer(wi), intent(in) :: nmax, mmax       ! m linear systems or size n
        real(wp), intent(in) :: rhs(nmax, 4)        ! RHS diagonals (#=5-1 because of normalization)
        real(wp), intent(in) :: u(mmax, nmax)       ! function u
        real(wp), intent(out) :: f(mmax, nmax)      ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1 ! rhs(1,1) contains 3. superdiagonal to allow for longer stencil at boundary
        f(:, n) = &
            +u(:, n) &
            + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4) + u(:, n + 3)*rhs(n, 1)

        n = 2
        f(:, n) = u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3)

        ! Interior points
        do n = 3, nmax - 2
            f(:, n) = u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                      + u(:, n) &
                      + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4)
        end do

        ! Boundary
        n = nmax - 1
        f(:, n) = u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3)

        n = nmax ! rhs(1,4) contains 3. subdiagonal to allow for longer stencil at boundary
        f(:, n) = u(:, n - 3)*rhs(n, 4) + u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                  + u(:, n)

        return
    end subroutine MatMul_5d

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. superdiagonal equal to 1
    subroutine MatMul_5d_antisym(nmax, mmax, r1, r2, r3, r4, r5, u, f, periodic, ibc)
        integer(wi), intent(in) :: nmax, mmax       ! m linear systems or size n
        real(wp), intent(in) :: r1(nmax), r2(nmax), r3(nmax), r4(nmax), r5(nmax)  ! RHS diagonals
        real(wp), intent(in) :: u(mmax, nmax)       ! function u
        real(wp), intent(out) :: f(mmax, nmax)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r5_loc     ! 2. off-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r5_loc = r5(3)      ! The first 2 equations, last 2 equations, are normalized differently

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nmax) &
                      + r5_loc*(u(:, 3) - u(:, nmax - 1))

            f(:, 2) = u(:, 3) - u(:, 1) &
                      + r5_loc*(u(:, 4) - u(:, nmax))

        else
            f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) &
                      + u(:, 4)*r1(1)   ! r1(1) contains 3. superdiagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 3, nmax - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nmax - 1) = u(:, nmax) - u(:, nmax - 2) &
                             + r5_loc*(u(:, 1) - u(:, nmax - 3))

            f(:, nmax) = u(:, 1) - u(:, nmax - 1) &
                         + r5_loc*(u(:, 2) - u(:, nmax - 2))

        else
            f(:, nmax - 1) = u(:, nmax - 3)*r1(nmax - 1) + u(:, nmax - 2)*r2(nmax - 1) + u(:, nmax - 1)*r3(nmax - 1) &
                             + u(:, nmax)*r4(nmax - 1)
            f(:, nmax) = u(:, nmax - 3)*r5(nmax) & ! r5(nmax) contains 3. subdiagonal to allow for longer stencil at boundary
                         + u(:, nmax - 2)*r1(nmax) + u(:, nmax - 1)*r2(nmax) + u(:, nmax)*r3(nmax)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nmax) = 0.0_wp

        end if

        return
    end subroutine MatMul_5d_antisym

    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. superdiagonal equal to 1
    subroutine MatMul_5d_sym(nmax, mmax, r0, r2, u, f)
        integer(wi), intent(in) :: nmax, mmax       ! m linear systems or size n
        real(wp), intent(in) :: r0, r2              ! diagonal and 2. off-diagonal
        real(wp), intent(in) :: u(mmax, nmax)       ! function u
        real(wp), intent(out) :: f(mmax, nmax)      ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1
        f(:, n) = r0*u(:, n) + u(:, n + 1) + u(:, nmax) &
                  + r2*(u(:, n + 2) + u(:, nmax - 1))

        n = 2
        f(:, n) = r0*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                  + r2*(u(:, n + 2) + u(:, nmax))

        ! Interior points
        do n = 3, nmax - 2
            f(:, n) = r0*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r2*(u(:, n + 2) + u(:, n - 2))
        end do

        ! Boundary
        n = nmax - 1
        f(:, n) = r0*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                  + r2*(u(:, 1) + u(:, n - 2))

        n = nmax
        f(:, n) = r0*u(:, n) + u(:, 1) + u(:, n - 1) &
                  + r2*(u(:, 2) + u(:, n - 2))

        return
    end subroutine MatMul_5d_sym

    ! Calculate f = B u, assuming B is antisymmetric hepta-diagonal with 1. superdiagonal equal to 1
    subroutine MatMul_7d_sym(nmax, mmax, r0, r2, r3, u, f)
        integer(wi), intent(in) :: nmax, mmax       ! m linear systems or size n
        real(wp), intent(in) :: r0, r2, r3          ! diagonal and 2. off-diagonal
        real(wp), intent(in) :: u(mmax, nmax)       ! function u
        real(wp), intent(out) :: f(mmax, nmax)      ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = r0*u(:, 1) + u(:, 2) + u(:, nmax) &
                  + r2*(u(:, 3) + u(:, nmax - 1)) &
                  + r3*(u(:, 4) + u(:, nmax - 2))

        f(:, 2) = r0*u(:, 2) + u(:, 3) + u(:, 1) &
                  + r2*(u(:, 4) + u(:, nmax)) &
                  + r3*(u(:, 5) + u(:, nmax - 1))

        f(:, 3) = r0*u(:, 3) + u(:, 4) + u(:, 2) &
                  + r2*(u(:, 5) + u(:, 1)) &
                  + r3*(u(:, 6) + u(:, nmax))

        ! Interior points
        do n = 4, nmax - 3
            f(:, n) = r0*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r2*(u(:, n + 2) + u(:, n - 2)) &
                      + r3*(u(:, n + 3) + u(:, n - 3))
        end do

        ! Boundary
        f(:, nmax - 2) = r0*u(:, nmax - 2) + u(:, nmax - 1) + u(:, nmax - 3) &
                         + r2*(u(:, nmax) + u(:, nmax - 4)) &
                         + r3*(u(:, 1) + u(:, nmax - 5))

        f(:, nmax - 1) = r0*u(:, nmax - 1) + u(:, nmax) + u(:, nmax - 2) &
                         + r2*(u(:, 1) + u(:, nmax - 3)) &
                         + r3*(u(:, 2) + u(:, nmax - 4))

        f(:, nmax) = r0*u(:, nmax) + u(:, 1) + u(:, nmax - 1) &
                     + r2*(u(:, 2) + u(:, nmax - 2)) &
                     + r3*(u(:, 3) + u(:, nmax - 3))

        return
    end subroutine MatMul_7d_sym

! #######################################################################
! #######################################################################
! TO BE IMPLEMENTED ACCORDING TO NOTES
    subroutine FDM_Bcs(lhs, ibc)
        real(wp), intent(inout) :: lhs(:, :)
        integer, intent(in) :: ibc

        integer(wi) id, nx

        id = size(lhs, 2)/2 + 1     ! central diagonal
        nx = size(lhs, 1)           ! # grid points

        if (any([BCS_ND, BCS_NN] == ibc)) then
            lhs(1, id:) = 0.0_wp
            lhs(1, id) = 1.0_wp
        end if

        if (any([BCS_DN, BCS_NN] == ibc)) then
            lhs(nx, :id) = 0.0_wp
            lhs(nx, id) = 1.0_wp
        end if

        return
    end subroutine FDM_Bcs

end module FDM_PROCS
