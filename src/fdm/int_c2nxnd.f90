!########################################################################
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme for non-uniform grids from Shukla and Zhong (2005).
!# in order to solve the BVP
!#
!#     u''_i - \lamba^2 u_i = f_i  N-2 eqns
!#     u_1 and u_N given           2   eqns
!#     Au'' = Bu                   N   eqns
!#
!# The system of N-2 eqns:
!#
!#                    (B - \lambda^2 A)u = Af = g
!#
!# is established in this routine, giving diagonals a-e and g (see notes).
!# See also routine fdm_c2n6nd.
!#
!# Solution array does not appear in this routine.
!#
!########################################################################

!########################################################################
!Left-hand side; pentadiagonal matrix of the linear system and f1 and f2
!########################################################################
subroutine INT_C2NXND_LHS_E(imax, x, ibc, lhs, rhs, lambda2, a, b, c, d, e, f1, f2)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(in) :: imax       ! original size; here using only 2:imax-1
    real(wp), intent(in) :: x(imax)
    integer, intent(in) :: ibc
    real(wp), intent(in) :: lhs(imax, 3)
    real(wp), intent(in) :: rhs(imax, 4)
    real(wp) lambda2
    real(wp), intent(out) :: a(imax), b(imax), c(imax), d(imax), e(imax)  ! diagonals
    real(wp), intent(out) :: f1(imax), f2(imax)      ! forcing term for the hyperbolic sine

! -------------------------------------------------------------------
    integer(wi) i
    real(wp) dummy1, dummy2, pprime, coef(4)
    integer(wi) stencil(4), size_of_stencil

! ###################################################################
! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system (array C22R)
! -------------------------------------------------------------------
    i = 2; dummy1 = lhs(2, 1)/lhs(1, 2)
    a(i) = 0.0_wp ! padding
    b(i) = 0.0_wp ! padding
    c(i) = 1.0_wp - lambda2*lhs(i, 2) - dummy1*(rhs(1, 3) - lambda2*lhs(1, 3))
    d(i) = rhs(i, 3) - lambda2*lhs(i, 3) - dummy1*rhs(1, 4)
    e(i) = rhs(i, 4) - dummy1*rhs(1, 1)

    i = 3
    a(i) = 0.0_wp ! padding
    b(i) = rhs(i, 2) - lambda2*lhs(i, 1)
    c(i) = 1.0_wp - lambda2*lhs(i, 2)
    d(i) = rhs(i, 3) - lambda2*lhs(i, 3)
    e(i) = rhs(i, 4)

    do i = 4, imax - 3
        a(i) = rhs(i, 1)
        b(i) = rhs(i, 2) - lambda2*lhs(i, 1)
        c(i) = 1.0_wp - lambda2*lhs(i, 2)
        d(i) = rhs(i, 3) - lambda2*lhs(i, 3)
        e(i) = rhs(i, 4)
    end do

    i = imax - 2
    a(i) = rhs(i, 1)
    b(i) = rhs(i, 2) - lambda2*lhs(i, 1)
    c(i) = 1.0_wp - lambda2*lhs(i, 2)
    d(i) = rhs(i, 3) - lambda2*lhs(i, 3)
    e(i) = 0.0_wp

    i = imax - 1; dummy2 = lhs(imax - 1, 3)/lhs(imax, 2)
    a(i) = rhs(i, 1) - dummy2*rhs(imax, 4)
    b(i) = rhs(i, 2) - lambda2*lhs(i, 1) - dummy2*rhs(imax, 1)
    c(i) = 1.0_wp - lambda2*lhs(i, 2) - dummy2*(rhs(imax, 2) - lambda2*lhs(imax, 1))
    d(i) = 0.0_wp ! padding
    e(i) = 0.0_wp ! padding

! -------------------------------------------------------------------
! Setting the RHS for the hyperbolic sine
! The minus sign in included here to save ops
! -------------------------------------------------------------------
    f1 = 0.0_wp ! b21R
    f1(1) = 1.0_wp    ! This element is simply the solution at imin of s(-)
    f1(2) = -(rhs(2, 2) - dummy1)
    f1(3) = -rhs(3, 1)

    f2 = 0.0_wp ! b2nR
    f2(imax - 2) = -rhs(imax - 2, 4)
    f2(imax - 1) = -(rhs(imax - 1, 3) - dummy2)
    f2(imax) = 1.0_wp ! This element is simply the solution at imax of s(+)

    ! -------------------------------------------------------------------
    ! Corrections to the BCS_DD to account for Neumann using third-order fdm for derivative at the boundary
    if (any([BCS_ND, BCS_NN] == ibc)) then
        ! Coefficients in FDM p'_1= c_1 p_1 + c_2 p_2 + c_3 p_3 + ...
        ! coef = [-1.5_wp, 2.0_wp, -0.5_wp, 0.0_wp]/dx(1,2)                       ! second-order, uniform
        ! coef = [-11.0_wp/6.0_wp, 3.0_wp, -1.5_wp, 1.0_wp/3.0_wp]/dx(1, 1)       ! third-order, uniform
        ! stencil(1:3) = [1, 2, 3]; size_of_stencil = 3                           ! second-order
        stencil(1:4) = [1, 2, 3, 4]; size_of_stencil = 4                        ! third-order
        coef(:) = 0.0_wp
        do i = 1, size_of_stencil
            coef(i) = Lag_p(x, 1, i, stencil(1:size_of_stencil))
        end do
        ! Data to calculate p_1 in terms of p_2, p_3, p_4 and p'_1
        c(1) = -coef(2)/coef(1)
        d(1) = -coef(3)/coef(1)
        e(1) = -coef(4)/coef(1)
        pprime = 1.0_wp/coef(1)
        ! Derived coefficients; see notes
        c(2) = c(2) - c(1)*f1(2)                                ! in reduced matrix; the minus sign comes from def of f1
        d(2) = d(2) - d(1)*f1(2)
        e(2) = e(2) - e(1)*f1(2)
        b(3) = b(3) - c(1)*f1(3)
        c(3) = c(3) - d(1)*f1(3)
        d(3) = d(3) - e(1)*f1(3)
        f1(:) = pprime*f1(:)                                    ! for particular solutions
    end if
    if (any([BCS_DN, BCS_NN] == ibc)) then
        ! Coefficients in FDM p'_n=c_1 p_n + c_2 p_{n-1} + c_3 p_{n-2} +...
        ! coef = [1.5_wp, -2.0_wp, 0.5_wp, 0.0_wp]/dx(imax, 1)                        ! second-order
        ! coef = [11.0_wp/6.0_wp, -3.0_wp, 1.5_wp, -1.0_wp/3.0_wp]/dx(imax, 1)        ! third-order
        ! stencil(1:3) = [imax, imax - 1, imax - 2]; size_of_stencil = 3              ! second-order
        stencil(1:4) = [imax, imax - 1, imax - 2, imax - 3]; size_of_stencil = 4    ! third-order
        coef(:) = 0.0_wp
        do i = 1, size_of_stencil
            coef(i) = Lag_p(x, imax, stencil(i), stencil(1:size_of_stencil))
        end do
        ! Data to calculate p_n in terms of p_{n-1}, p_{n-2} and p'_n
        c(imax) = -coef(2)/coef(1)
        b(imax) = -coef(3)/coef(1)
        a(imax) = -coef(4)/coef(1)
        pprime = 1.0_wp/coef(1)
        ! Derived coefficients; see notes
        a(imax - 1) = a(imax - 1) - a(imax)*f2(imax - 1)        ! in reduced matrix; the minus sign comes from def of f2
        b(imax - 1) = b(imax - 1) - b(imax)*f2(imax - 1)
        c(imax - 1) = c(imax - 1) - c(imax)*f2(imax - 1)
        b(imax - 2) = b(imax - 2) - a(imax)*f2(imax - 2)
        c(imax - 2) = c(imax - 2) - b(imax)*f2(imax - 2)
        d(imax - 2) = d(imax - 2) - c(imax)*f2(imax - 2)
        f2(:) = pprime*f2(:)                                    ! for particular solutions
    end if

    return
contains
    function Lag_p(x, j, i, idx) result(f)                     ! 1. derivative of Lagrange polynomials on idx(:) around i
        use TLAB_CONSTANTS
        implicit none

        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f, num, den, dummy
        integer(wi) k, m

        den = 1.0_wp
        num = 0.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                dummy = 1.0_wp
                do m = 1, size(idx)
                    if (idx(m) /= i .and. m /= k) then
                        dummy = dummy*(x(j) - x(idx(m)))
                    end if
                end do
                num = num + dummy
                den = den*(x(i) - x(idx(k)))
            end if
        end do
        f = num/den

        return
    end function

end subroutine INT_C2NXND_LHS_E

! #######################################################################
! Right-hand side; mmax forcing terms at the same time
! #######################################################################
subroutine INT_C2NXND_RHS(imax, mmax, lhs, f, g)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(in) :: imax, mmax
    real(wp), intent(in) :: lhs(imax, 3)
    real(wp), intent(in) :: f(mmax, imax)
    real(wp), intent(out) :: g(mmax, imax)

! -------------------------------------------------------------------
    integer(wi) i
    real(wp) dummy

! ###################################################################
    g(:, 1) = 0.0_wp ! This element is simply the solution at imin of p(0)

    i = 2
    dummy = lhs(i, 2) - lhs(2, 1)/lhs(1, 2)*lhs(1, 3)
    g(:, i) = f(:, i)*dummy + f(:, i + 1)*lhs(i, 3)

    do i = 3, imax - 2 ! Interior points
        g(:, i) = f(:, i - 1)*lhs(i, 1) + f(:, i)*lhs(i, 2) + f(:, i + 1)*lhs(i, 3)
    end do

    i = imax - 1
    dummy = lhs(i, 2) - lhs(imax - 1, 3)/lhs(imax, 2)*lhs(imax, 1)
    g(:, i) = f(:, i - 1)*lhs(i, 1) + f(:, i)*dummy

    g(:, imax) = 0.0_wp ! This element is simply the solution at imax of p(0)

    return
end subroutine INT_C2NXND_RHS
