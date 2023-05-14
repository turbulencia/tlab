#define C_12D11_L .1090909090909090e+1_wp
#define C_51D22_L .2318181818181818e+1_wp
#define C_03D44_L .6818181818181818e-1_wp
#define C_02D11_L .1818181818181818e+0_wp

!########################################################################
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic,
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
!# Interior points 6th-order according to Eq. 2.2.7.
!# The second point from Eq. 2.2.6 forth-order (b=0).
!# The first point from third-order biased Eq. 4.3.1.
!#
!# Solution array does not appear in this routine.
!#
!########################################################################

!########################################################################
!Left-hand side; pentadiagonal matrix of the linear system and f1 and f2
!########################################################################
subroutine INT_C2N6_LHS_E(imax, dx, ibc, lambda2, a, b, c, d, e, f1, f2)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(in) :: imax       ! original size; here using only 2:imax-1
    real(wp), intent(in) :: dx(imax, 2)
    integer, intent(in) :: ibc
    real(wp), intent(in) :: lambda2
    real(wp), intent(out) :: a(imax), b(imax), c(imax), d(imax), e(imax)  ! diagonals
    real(wp), intent(out) :: f1(imax), f2(imax)     ! forcing term for the hyperbolic sine

! -------------------------------------------------------------------
    integer(wi) i

! ###################################################################
! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
    a(2) = 0.0_wp ! padding
    b(2) = 0.0_wp ! padding
    c(2) = -3.0_wp - lambda2*dx(2, 1)*dx(2, 1)
    d(2) = 3.0_wp + lambda2*dx(3, 1)*dx(3, 1)
    e(2) = -1.0_wp
! fourth-order centered
    a(3) = 0.0_wp ! padding
    b(3) = C_12D11_L - C_02D11_L*lambda2*dx(2, 1)*dx(2, 1)
    c(3) = -C_51D22_L - lambda2*dx(3, 1)*dx(3, 1)
    d(3) = C_12D11_L - C_02D11_L*lambda2*dx(4, 1)*dx(4, 1)
    e(3) = C_03D44_L
! sixth-order centered
    do i = 4, imax - 3
        a(i) = C_03D44_L
        b(i) = C_12D11_L - C_02D11_L*lambda2*dx(i - 1, 1)*dx(i - 1, 1)
        c(i) = -C_51D22_L - lambda2*dx(i, 1)*dx(i, 1)
        d(i) = C_12D11_L - C_02D11_L*lambda2*dx(i + 1, 1)*dx(i + 1, 1)
        e(i) = C_03D44_L
    end do
! fourth-order centered
    a(imax - 2) = C_03D44_L
    b(imax - 2) = C_12D11_L - C_02D11_L*lambda2*dx(imax - 3, 1)*dx(imax - 3, 1)
    c(imax - 2) = -C_51D22_L - lambda2*dx(imax - 2, 1)*dx(imax - 2, 1)
    d(imax - 2) = C_12D11_L - C_02D11_L*lambda2*dx(imax - 1, 1)*dx(imax - 1, 1)
    e(imax - 2) = 0.0_wp ! padding
! third-order biased
    a(imax - 1) = -1.0_wp
    b(imax - 1) = 3.0_wp + lambda2*dx(imax - 2, 1)*dx(imax - 2, 1)
    c(imax - 1) = -3.0_wp - lambda2*dx(imax - 1, 1)*dx(imax - 1, 1)
    d(imax - 1) = 0.0_wp ! padding
    e(imax - 1) = 0.0_wp ! padding

! -------------------------------------------------------------------
! Setting the RHS for the hyperbolic sine
! The minus sign in included here to save ops
! -------------------------------------------------------------------
    f1 = 0.0_wp
    f1(1) = 1.0_wp ! This element is simply the solution at imin of s(-)
    f1(2) = -1.0_wp
    f1(3) = -C_03D44_L

    f2 = 0.0_wp
    f2(imax - 2) = -C_03D44_L
    f2(imax - 1) = -1.0_wp
    f2(imax) = 1.0_wp ! This element is simply the solution at imax of s(+)

    ! -------------------------------------------------------------------
!   Corrections to the BCS_DD to account for Neumann using second-order fdm for derivative at the boundary
    if (any([BCS_ND, BCS_NN] == ibc)) then
        c(1) = 4.0_wp/3.0_wp                                ! To calculate p_1 in terms of p_2, p_3
        d(1) = -1.0_wp/3.0_wp
        e(1) = 0.0_wp
        c(2) = c(2) + 4.0_wp/3.0_wp                         ! in reduced matrix
        d(2) = d(2) - 1.0_wp/3.0_wp
        b(3) = b(3) + 1.0_wp/11.0_wp
        c(3) = c(3) - 1.0_wp/44.0_wp
        f1(:) = -2.0_wp/3.0_wp*dx(1, 1)*f1(:)               ! for particular solutions
    end if
    if (any([BCS_DN, BCS_NN] == ibc)) then
        c(imax - 2) = c(imax - 2) - 1.0_wp/44.0_wp          ! in reduced matrix
        d(imax - 2) = d(imax - 2) + 1.0_wp/11.0_wp
        b(imax - 1) = b(imax - 1) - 1.0_wp/3.0_wp
        c(imax - 1) = c(imax - 1) + 4.0_wp/3.0_wp
        f2(:) = 2.0_wp/3.0_wp*dx(imax, 1)*f2(:)             ! for particular solutions
        a(imax) = 0.0_wp
        b(imax) = -1.0_wp/3.0_wp                            ! To calculate p_n in terms of p_{n-1}, p_{n-2}
        c(imax) = 4.0_wp/3.0_wp
    end if

    return
end subroutine INT_C2N6_LHS_E

! #######################################################################
! Right-hand side; jkmax forcing terms at the same time
! #######################################################################
subroutine INT_C2N6_RHS(imax, jkmax, dx, f, g)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(in) :: imax, jkmax
    real(wp), intent(in) :: dx(imax, 2)
    real(wp), intent(in) :: f(jkmax, imax)
    real(wp), intent(out) :: g(jkmax, imax)

! -------------------------------------------------------------------
    integer(wi) i

! ###################################################################
! Boundary conditions
    g(:, 1) = 0.0_wp ! This element is simply the solution at imin of p(0)
    g(:, 2) = f(:, 2)*dx(2, 1)*dx(2, 1) - f(:, 3)*dx(3, 1)*dx(3, 1)
    g(:, imax - 1) = f(:, imax - 1)*dx(imax - 1, 1)*dx(imax - 1, 1) - f(:, imax - 2)*dx(imax - 2, 1)*dx(imax - 2, 1)
    g(:, imax) = 0.0_wp ! This element is simply the solution at imax of p(0)

! Interior points
    do i = 3, imax - 2
        g(:, i) = f(:, i)*dx(i, 1)*dx(i, 1) + &
                  C_02D11_L*(f(:, i - 1)*dx(i - 1, 1)*dx(i - 1, 1) + f(:, i + 1)*dx(i + 1, 1)*dx(i + 1, 1))
    end do

    return
end subroutine INT_C2N6_RHS
