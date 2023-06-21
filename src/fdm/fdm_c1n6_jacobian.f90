!########################################################################
!# Compact FDMs for non-uniform grids from Lele, JCP, 1992, usign Jacobian
!#
!# f'_i + a_1(f'_{i-1}+f'_{i+1}) + a_2(f'_{i-1}+f'_{i+1})
!#   = b_1(f_{i+1}-f_{i-1}) + b_2(f_{i+2}-f_{i-2}) + b_3(f_{i+3}-f_{i-3})
!#
!# Here we assume b_3 = 0
!#
!# System normalized s.t. 1. off-diagonal in B is 1.
!#
!########################################################################
subroutine FDM_C1N6_Jacobian(nmax, dx, lhs, rhs, coef, periodic)
    use TLAB_CONSTANTS, only: wp, wi
    implicit none

    integer(wi), intent(in) :: nmax
    real(wp), intent(in) :: dx(nmax)
    real(wp), intent(out) :: lhs(nmax, 3)   ! LHS diagonals
    real(wp), intent(out) :: rhs(nmax, 5)   ! RHS diagonals; assuming b_3 = 0
    real(wp), intent(out) :: coef(5)        ! a_1, a_2, b_1, b_2, b_3
    logical, intent(in) :: periodic

    integer(wi) n

    ! #######################################################################
    ! Interior points according to Eq. 2.1.7 (\alpha=1/3), 6th order approximation
    coef(1:2) = [1.0_wp/3.0_wp, 0.0_wp]                     ! a_1, a_2
    coef(3:5) = [7.0_wp/9.0_wp, 1.0_wp/36.0_wp, 0.0_wp]     ! b_1, b_2, b_3

    ! lhs diagonals
    lhs(:, 2) = 1.0_wp              ! diagonal
    lhs(:, [1, 3]) = coef(1)        ! 1. off-diagonal

    ! rhs diagonals
    rhs(:, 1) = -coef(4)            ! 2. subdiagonal
    rhs(:, 2) = -coef(3)            ! 1. subdiagonal
    rhs(:, 3) = 0.0_wp              ! diagonal
    rhs(:, 4) = coef(3)             ! 1. superdiagonal
    rhs(:, 5) = coef(4)             ! 2. superdiagonal

    ! boundaries
    if (.not. periodic) then
        ! 3. order in Eq. 4.1.3 with \alpha=2
        n = 1
        lhs(n, :) = [0.0_wp, 1.0_wp, 2.0_wp]
        rhs(n, :) = [0.0_wp, 0.0_wp, -15.0_wp/6.0_wp, 2.0_wp, 0.5_wp]

        ! 5. order in Carpenter, Gottlieb and Aberbanel, JCP, 1993, Eq. 93
        ! who study the effect of boundary points on stability
        n = 2
        lhs(n, :) = [1.0_wp/6.0_wp, 1.0_wp, 0.5_wp]
        rhs(n, :) = [0.0_wp, -5.0_wp/9.0_wp, -0.5_wp, 1.0_wp, 1.0_wp/18.0_wp]
        ! 4. order in Eq. 2.1.6 with \alpha=1/4.
        ! lhs(n, :) = [0.25_wp, 1.0_wp, 0.25_wp]
        ! rhs(n, :) = [0.0_wp, -0.75_wp, 0.0_wp, 0.75_wp, 0.0_wp]

        ! symmetry propertie to define values at end
        n = nmax - 1
        lhs(n, :) = lhs(2, 3:1:-1)
        rhs(n, :) = -rhs(2, 5:1:-1)

        n = nmax
        lhs(n, :) = lhs(1, 3:1:-1)
        rhs(n, :) = -rhs(1, 5:1:-1)
    end if

    ! multiply by the Jacobian
    lhs(:, 1) = lhs(:, 1)*cshift(dx(:), -1)
    lhs(:, 2) = lhs(:, 2)*dx(:)
    lhs(:, 3) = lhs(:, 3)*cshift(dx(:), +1)

    ! normalize such the coefficent in 1. off-diagonal of rhs is 1
    lhs(:, :) = lhs(:, :)/coef(3)
    rhs(:, :) = rhs(:, :)/coef(3)

    return
end subroutine FDM_C1N6_Jacobian
