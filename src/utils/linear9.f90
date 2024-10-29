!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/08/07 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Solve LU of nonadiagonal system
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
subroutine NONADFS(nmax, inorm, a, b, c, d, e, f, g, h, i)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi) nmax, inorm
    real(wp), dimension(nmax), intent(INOUT) :: a, b, c, d, e, f, g, h, i

! -----------------------------------------------------------------------
    integer(wi) n

! #######################################################################
    a(1) = 1.0_wp ! padding
    b(1) = 1.0_wp ! padding
    c(1) = 1.0_wp ! padding
    if (inorm == 1) then
        d(1) = 1.0_wp/e(1) ! padding, and used in nonadss to normalize 1st eqn.
        e(1) = 1.0_wp
        f(1) = f(1)*d(1)
        g(1) = g(1)*d(1)
        h(1) = h(1)*d(1)
        i(1) = i(1)*d(1)
    else
        d(1) = 1.0_wp ! padding
    end if

!  n = 2
    a(2) = 1.0_wp ! padding
    b(2) = 1.0_wp ! padding
    c(2) = 1.0_wp ! padding
    d(2) = d(2)/e(1)
    e(2) = e(2) - d(2)*f(1)
    f(2) = f(2) - d(2)*g(1)
    g(2) = g(2) - d(2)*h(1)
    h(2) = h(2) - d(2)*i(1)

    n = 3
    a(n) = 1.0_wp ! padding
    b(n) = 1.0_wp ! padding
    c(n) = c(n)/e(n - 2)
    d(n) = (d(n) - c(n)*f(n - 2))/e(n - 1)
    e(n) = e(n) - d(n)*f(n - 1) - c(n)*g(n - 2)
    f(n) = f(n) - d(n)*g(n - 1) - c(n)*h(n - 2)
    g(n) = g(n) - d(n)*h(n - 1) - c(n)*i(n - 2)
    h(n) = h(n) - d(n)*i(n - 1)

    n = 4
    a(n) = 1.0_wp ! padding
    b(n) = b(n)/e(n - 3)
    c(n) = (c(n) - b(n)*f(n - 3))/e(n - 2)
    d(n) = (d(n) - c(n)*f(n - 2) - b(n)*g(n - 3))/e(n - 1)
    e(n) = e(n) - d(n)*f(n - 1) - c(n)*g(n - 2) - b(n)*h(n - 3)
    f(n) = f(n) - d(n)*g(n - 1) - c(n)*h(n - 2) - b(n)*i(n - 3)
    g(n) = g(n) - d(n)*h(n - 1) - c(n)*i(n - 2)
    h(n) = h(n) - d(n)*i(n - 1)

! -----------------------------------------------------------------------
    do n = 5, nmax - 3
        a(n) = a(n)/e(n - 4)
        b(n) = (b(n) - a(n)*f(n - 4))/e(n - 3)
        c(n) = (c(n) - b(n)*f(n - 3) - a(n)*g(n - 4))/e(n - 2)
        d(n) = (d(n) - c(n)*f(n - 2) - b(n)*g(n - 3) - a(n)*h(n - 4))/e(n - 1)
        e(n) = e(n) - d(n)*f(n - 1) - c(n)*g(n - 2) - b(n)*h(n - 3) - a(n)*i(n - 4)
        f(n) = f(n) - d(n)*g(n - 1) - c(n)*h(n - 2) - b(n)*i(n - 3)
        g(n) = g(n) - d(n)*h(n - 1) - c(n)*i(n - 2)
        h(n) = h(n) - d(n)*i(n - 1)
    end do
    i(n - 1) = 1.0_wp ! padding

! -----------------------------------------------------------------------
    n = nmax - 2
    a(n) = a(n)/e(n - 4)
    b(n) = (b(n) - a(n)*f(n - 4))/e(n - 3)
    c(n) = (c(n) - b(n)*f(n - 3) - a(n)*g(n - 4))/e(n - 2)
    d(n) = (d(n) - c(n)*f(n - 2) - b(n)*g(n - 3) - a(n)*h(n - 4))/e(n - 1)
    e(n) = e(n) - d(n)*f(n - 1) - c(n)*g(n - 2) - b(n)*h(n - 3) - a(n)*i(n - 4)
    f(n) = f(n) - d(n)*g(n - 1) - c(n)*h(n - 2) - b(n)*i(n - 3)
    g(n) = g(n) - d(n)*h(n - 1) - c(n)*i(n - 2)
    h(n) = 1.0_wp ! padding
    i(n) = 1.0_wp ! padding

    n = nmax - 1
    a(n) = a(n)/e(n - 4)
    b(n) = (b(n) - a(n)*f(n - 4))/e(n - 3)
    c(n) = (c(n) - b(n)*f(n - 3) - a(n)*g(n - 4))/e(n - 2)
    d(n) = (d(n) - c(n)*f(n - 2) - b(n)*g(n - 3) - a(n)*h(n - 4))/e(n - 1)
    e(n) = e(n) - d(n)*f(n - 1) - c(n)*g(n - 2) - b(n)*h(n - 3) - a(n)*i(n - 4)
    f(n) = f(n) - d(n)*g(n - 1) - c(n)*h(n - 2) - b(n)*i(n - 3)
    g(n) = 1.0_wp ! padding
    h(n) = 1.0_wp ! padding
    i(n) = 1.0_wp ! padding

    n = nmax
    a(n) = a(n)/e(n - 4)
    b(n) = (b(n) - a(n)*f(n - 4))/e(n - 3)
    c(n) = (c(n) - b(n)*f(n - 3) - a(n)*g(n - 4))/e(n - 2)
    d(n) = (d(n) - c(n)*f(n - 2) - b(n)*g(n - 3) - a(n)*h(n - 4))/e(n - 1)
    e(n) = e(n) - d(n)*f(n - 1) - c(n)*g(n - 2) - b(n)*h(n - 3) - a(n)*i(n - 4)
    f(n) = 1.0_wp ! padding
    g(n) = 1.0_wp ! padding
    h(n) = 1.0_wp ! padding
    i(n) = 1.0_wp ! padding

    return
end subroutine NONADFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
subroutine NONADSS(nmax, len, a, b, c, d, e, f, g, h, i, frc)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi) nmax, len
    real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e, f, g, h, i
    real(wp), dimension(len, nmax), intent(INOUT) :: frc

! -----------------------------------------------------------------------
    integer(wi) n, ij

! #######################################################################
! -----------------------------------------------------------------------
! Solve Ly=frc, forward
! -----------------------------------------------------------------------
    do ij = 1, len
        frc(ij, 1) = frc(ij, 1)*d(1) ! Normalize first eqn. See NONADFS
        frc(ij, 2) = frc(ij, 2) - frc(ij, 1)*d(2)
        frc(ij, 3) = frc(ij, 3) - frc(ij, 2)*d(3) - frc(ij, 1)*c(3)
        frc(ij, 4) = frc(ij, 4) - frc(ij, 3)*d(4) - frc(ij, 2)*c(4) - frc(ij, 1)*b(4)
    end do

    do n = 5, nmax
        do ij = 1, len
            frc(ij, n) = frc(ij, n) - frc(ij, n - 1)*d(n) - frc(ij, n - 2)*c(n) - frc(ij, n - 3)*b(n) - frc(ij, n - 4)*a(n)
        end do
    end do

! -----------------------------------------------------------------------
! Solve Ux=y, backward
! -----------------------------------------------------------------------
    do ij = 1, len
        frc(ij, nmax) = frc(ij, nmax) &
                        /e(nmax)
        frc(ij, nmax - 1) = (frc(ij, nmax - 1) - frc(ij, nmax)*f(nmax - 1) &
                             )/e(nmax - 1)
        frc(ij, nmax - 2) = (frc(ij, nmax - 2) - frc(ij, nmax - 1)*f(nmax - 2) - frc(ij, nmax)*g(nmax - 2) &
                             )/e(nmax - 2)
        frc(ij, nmax - 3) = (frc(ij, nmax - 3) - frc(ij, nmax - 2)*f(nmax - 3) - frc(ij, nmax - 1)*g(nmax - 3) - frc(ij, nmax)*h(nmax - 3) &
                             )/e(nmax - 3)
    end do

    do n = nmax - 4, 1, -1
        do ij = 1, len
            frc(ij, n) = (frc(ij, n) - frc(ij, n + 1)*f(n) - frc(ij, n + 2)*g(n) - frc(ij, n + 3)*h(n) - frc(ij, n + 4)*i(n))/e(n)
        end do
    end do

    return
end subroutine NONADSS
