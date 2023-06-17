!########################################################################
!#
!# Initialize the solver for the BVP
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
!########################################################################

!########################################################################
!Left-hand side; pentadiagonal matrix of the linear system and f1 and f2
!########################################################################
subroutine INT_C2NX_LHS_E(imax, x, ibc, lhs, rhs, lambda2, a, b, c, d, e, f1, f2)
    use TLAB_CONSTANTS
    use FDM_PROCS
    implicit none

    integer(wi), intent(in) :: imax         ! original size; here using only 2:imax-1
    real(wp), intent(in) :: x(imax)         ! grid nodes
    integer, intent(in) :: ibc              ! Boundary condition, BCS_DD, BCS_DN, BCS_ND, BCS_NN
    real(wp), intent(in) :: lhs(imax, 3)    ! matrix A, tridiagonal
    real(wp), intent(in) :: rhs(imax, 4)    ! matrix B, with unitary central diagonal
    real(wp) lambda2                        ! system constatn
    real(wp), intent(out) :: a(imax), b(imax), c(imax), d(imax), e(imax)  ! diagonals
    real(wp), intent(out) :: f1(imax), f2(imax)      ! forcing term for the hyperbolic sine

! -------------------------------------------------------------------
    integer(wi) i
    real(wp) dummy1, dummy2, pprime, coef(5)

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
        ! Coefficients in FDM p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
        coef(:) = 0.0_wp
        ! coef(1:3) = coef_e1n2_biased(x, 1)                ! second-order
        coef(1:4) = coef_e1n3_biased(x, 1)                ! third-order
        ! coef(1:5) = coef_c1n4_biased(x, 1)                ! fourth-order

        ! Data to calculate p_1 in terms of p_2, p_3, p_4 and p'_1
        c(1) = -(coef(2) + lambda2*coef(5))/coef(1)         ! d + lambda^2h^2 e in notes, e only 1 component
        d(1) = -coef(3)/coef(1)
        e(1) = -coef(4)/coef(1)
        pprime = 1.0_wp/coef(1)
        b(1) = -coef(5)/coef(1)                             ! coefficient e for p''_2

        ! Derived coefficients; see notes
        c(2) = c(2) - c(1)*f1(2)                            ! in reduced C matrix; the minus sign comes from def of f1
        d(2) = d(2) - d(1)*f1(2)
        e(2) = e(2) - e(1)*f1(2)
        b(3) = b(3) - c(1)*f1(3)
        c(3) = c(3) - d(1)*f1(3)
        d(3) = d(3) - e(1)*f1(3)

        b(2) = b(1)*f1(2)                                   ! in reduced A matrix; the plus sign comes from def of f1
        a(3) = b(1)*f1(3)

        f1(:) = pprime*f1(:)                                ! for particular solutions

    end if
    if (any([BCS_DN, BCS_NN] == ibc)) then
        ! Coefficients in FDM p'_n = b_1 p_n + b_2 p_{n-1} + b_3 p_{n-2} +...
        coef(:) = 0.0_wp
        ! coef(1:3) = coef_e1n2_biased(x, imax, backwards=.true.)
        coef(1:4) = coef_e1n3_biased(x, imax, backwards=.true.)
        ! coef(1:5) = coef_c1n4_biased(x, imax, backwards=.true.)

        ! Data to calculate p_n in terms of p_{n-1}, p_{n-2} and p'_n
        c(imax) = -(coef(2) + lambda2*coef(5))/coef(1)
        b(imax) = -coef(3)/coef(1)
        a(imax) = -coef(4)/coef(1)
        pprime = 1.0_wp/coef(1)
        e(imax) = -coef(5)/coef(1)                          ! coefficient for p''_{n-1}

        ! Derived coefficients; see notes
        a(imax - 1) = a(imax - 1) - a(imax)*f2(imax - 1)    ! in reduced C matrix; the minus sign comes from def of f2
        b(imax - 1) = b(imax - 1) - b(imax)*f2(imax - 1)
        c(imax - 1) = c(imax - 1) - c(imax)*f2(imax - 1)
        b(imax - 2) = b(imax - 2) - a(imax)*f2(imax - 2)
        c(imax - 2) = c(imax - 2) - b(imax)*f2(imax - 2)
        d(imax - 2) = d(imax - 2) - c(imax)*f2(imax - 2)

        d(imax - 1) = e(imax)*f2(imax - 1)                  ! in reduced A matrix; the plus sign comes from def of f2
        e(imax - 2) = e(imax)*f2(imax - 2)

        f2(:) = pprime*f2(:)                                ! for particular solutions

    end if

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

end subroutine INT_C2NX_LHS_E

! #######################################################################
! Right-hand side; mmax forcing terms at the same time
! #######################################################################
subroutine INT_C2NX_RHS(imax, mmax, lhs, f, g)
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
end subroutine INT_C2NX_RHS
