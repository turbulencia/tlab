!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/05/21 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Solve LU of heptadiagonal system
!#
!# First equation is divided by d1 to get 1 in the first diagonal elem.
!#
!########################################################################
!# ARGUMENTS
!#
!# nmax    In     Size of the pentadiagonal system
!# len     In     Number of simultaneous systems to be solved
!# frc     In     Array with forcing term
!#         Out    Array with solution
!#
!########################################################################

! #######################################################################
! LU factorization stage
! #######################################################################
subroutine HEPTADFS(nmax, a, b, c, d, e, f, g)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi) nmax
    real(wp), dimension(nmax), intent(INOUT) :: a, b, c, d, e, f, g

! -----------------------------------------------------------------------
    integer(wi) n

! #######################################################################
    g(1) = g(1)/d(1)
    f(1) = f(1)/d(1)
    e(1) = e(1)/d(1)
    c(1) = 1.0_wp/d(1) ! padding, and used in heptadss to normalize 1st eqn.
!   b(1) = 1.0_wp      ! padding
!   a(1) = 1.0_wp      ! padding
    d(1) = 1.0_wp

!   a(2) = 1.0_wp ! padding
!   b(2) = 1.0_wp ! padding
    c(2) = c(2)/d(1)
    d(2) = d(2) - c(2)*e(1)
    e(2) = e(2) - c(2)*f(1)
    f(2) = f(2) - c(2)*g(1)

!   a(3) = 1.0_wp ! padding
    b(3) = b(3)/d(1)
    c(3) = (c(3) - b(3)*e(1))/d(2)
    d(3) = d(3) - c(3)*e(2) - b(3)*f(1)
    e(3) = e(3) - c(3)*f(2) - b(3)*g(1)
    f(3) = f(3) - c(3)*g(2)

    do n = 4, nmax - 2
        a(n) = a(n)/d(n - 3)
        b(n) = (b(n) - a(n)*e(n - 3))/d(n - 2)
        c(n) = (c(n) - b(n)*e(n - 2) - a(n)*f(n - 3))/d(n - 1)
        d(n) = d(n) - c(n)*e(n - 1) - b(n)*f(n - 2) - a(n)*g(n - 3)
        e(n) = e(n) - c(n)*f(n - 1) - b(n)*g(n - 2)
        f(n) = f(n) - c(n)*g(n - 1)
    end do
!   g(n-1) = 1.0_wp ! padding

    n = nmax - 1
    a(n) = a(n)/d(n - 3)
    b(n) = (b(n) - a(n)*e(n - 3))/d(n - 2)
    c(n) = (c(n) - b(n)*e(n - 2) - a(n)*f(n - 3))/d(n - 1)
    d(n) = d(n) - c(n)*e(n - 1) - b(n)*f(n - 2) - a(n)*g(n - 3)
    e(n) = e(n) - c(n)*f(n - 1) - b(n)*g(n - 2)
!   f(n) = 1.0_wp ! padding
!   g(n) = 1.0_wp ! padding

    n = nmax
    a(n) = a(n)/d(n - 3)
    b(n) = (b(n) - a(n)*e(n - 3))/d(n - 2)
    c(n) = (c(n) - b(n)*e(n - 2) - a(n)*f(n - 3))/d(n - 1)
    d(n) = d(n) - c(n)*e(n - 1) - b(n)*f(n - 2) - a(n)*g(n - 3)
!   e(n) = 1.0_wp ! padding
!   f(n) = 1.0_wp ! padding
!   g(n) = 1.0_wp ! padding

    return
end subroutine HEPTADFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
subroutine HEPTADSS(nmax, len, a, b, c, d, e, f, g, frc)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi) nmax, len
    real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e, f, g
    real(wp), dimension(len, nmax), intent(INOUT) :: frc

! -----------------------------------------------------------------------
    integer(wi) n, ij

! #######################################################################
! -----------------------------------------------------------------------
! Solve Ly=frc, forward
! -----------------------------------------------------------------------
    do ij = 1, len
        frc(ij, 1) = frc(ij, 1)*c(1) ! Normalize first eqn. See HEPTADFS
        frc(ij, 2) = frc(ij, 2) - frc(ij, 1)*c(2)
        frc(ij, 3) = frc(ij, 3) - frc(ij, 2)*c(3) - frc(ij, 1)*b(3)
    end do

    do n = 4, nmax
        do ij = 1, len
            frc(ij, n) = frc(ij, n) - frc(ij, n - 1)*c(n) - frc(ij, n - 2)*b(n) - frc(ij, n - 3)*a(n)
        end do
    end do

! -----------------------------------------------------------------------
! Solve Ux=y, backward
! -----------------------------------------------------------------------
    do ij = 1, len
        frc(ij, nmax) = frc(ij, nmax)/d(nmax)
        frc(ij, nmax - 1) = (frc(ij, nmax - 1) - frc(ij, nmax)*e(nmax - 1))/d(nmax - 1)
        frc(ij, nmax - 2) = (frc(ij, nmax - 2) - frc(ij, nmax - 1)*e(nmax - 2) - frc(ij, nmax)*f(nmax - 2))/d(nmax - 2)
    end do

    do n = nmax - 3, 1, -1
        do ij = 1, len
            frc(ij, n) = (frc(ij, n) - frc(ij, n + 1)*e(n) - frc(ij, n + 2)*f(n) - frc(ij, n + 3)*g(n))/d(n)
        end do
    end do

    return
end subroutine HEPTADSS
