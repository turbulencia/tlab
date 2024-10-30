!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2010/12/20 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 6th order pentadiagonal compact scheme by JCP Lele 1992, non-periodic.
!# Interior points according to Eq. 2.1.10 and closer to the boundaries
!# with a tridiagonal scheme (alpha=1/3) according to Eq. 2.1.7.
!#
!# The following IVP is solved
!#
!#     u'_i + \lamba u_i = h_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# The system of N-1 eqns:
!#
!#                    (B + \lambda A)u = Ah = l
!#
!# is established in this routine, giving diagonals a-g and l (see notes).
!#
!# Carpenter, Gottlieb and Aberbanel, JCP, 1993, study the effect of
!# boundary points on stability. Scheme 3-5-2xtri6c--penta6c--2xtri6c-5-3
!# is implmented.
!#
!# Solution array does not appear in this routine.
!#
!########################################################################
!# ARGUMENTS
!#
!# jkmax       In    number of systems to solve
!# ibc         In    BCs: 1 u_1 given
!#                        2 u_N given
!#
!# h           Out   forcing term for the exponential
!#
!########################################################################

!########################################################################
!Left-hand side; heptadiagonal matrix of the linear system and h
!########################################################################
subroutine INT_C1N6M_LHS_E(imax, ibc, dx, lambda, a, b, c, d, e, f, g, h)
    use TLab_Constants, only: wp, wi
    use FDM_PROCS, only: C1N6M_ALPHA, C1N6M_BETA
    use FDM_PROCS, only: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6

    implicit none

    real(wp) lambda
    integer(wi), intent(IN) :: imax, ibc
    real(wp), dimension(imax), intent(IN) :: dx
    real(wp), dimension(imax), intent(OUT) :: a, b, c, d, e, f, g, h

! -------------------------------------------------------------------
    integer(wi) :: i
    real(wp) :: c0136, c1418, c0103, c0104

! ###################################################################
    c0136 = 1.0_wp/36.0_wp
    c1418 = 14.0_wp/18.0_wp
    c0103 = 1.0_wp/3.0_wp
    c0104 = 1.0_wp/4.0_wp

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
    a(1) = 0.0_wp
    b(1) = 0.0_wp
    c(1) = 0.0_wp
    d(1) = -5.0_wp/2.0_wp + lambda*dx(1)
    e(1) = 2.0_wp + lambda*2.0_wp*dx(2)
    f(1) = 0.5_wp
    g(1) = 0.0_wp
! fifth-order biased
    a(2) = 0.0_wp
    b(2) = 0.0_wp
    c(2) = -5.0_wp/9.0_wp + lambda*1.0_wp/6.0_wp*dx(1)
    d(2) = -1.0_wp/2.0_wp + lambda*dx(2)
    e(2) = 1.0_wp + lambda*1.0_wp/2.0_wp*dx(3)
    f(2) = 1.0_wp/18.0_wp
    g(2) = 0.0_wp
! sixth-order centered (alpha=1/3)
    do i = 3, 4 ! DO i = 3,3
        a(i) = 0.0_wp
        b(i) = -c0136
        c(i) = -c1418 + lambda*c0103*dx(i - 1)
        d(i) = 0.0_wp + lambda*dx(i)
        e(i) = c1418 + lambda*c0103*dx(i + 1)
        f(i) = c0136
        g(i) = 0.0_wp
    end do
! sixth-order modified centered
    do i = 5, imax - 4 ! DO i = 4,imax-3
        a(i) = -C1N6M_CD6
        b(i) = -C1N6M_BD4 + lambda*C1N6M_BETA*dx(i - 2)
        c(i) = -C1N6M_AD2 + lambda*C1N6M_ALPHA*dx(i - 1)
        d(i) = 0.0_wp + lambda*dx(i)
        e(i) = C1N6M_AD2 + lambda*C1N6M_ALPHA*dx(i + 1)
        f(i) = C1N6M_BD4 + lambda*C1N6M_BETA*dx(i + 2)
        g(i) = C1N6M_CD6
    end do
! sixth-order centered (alpha=1/3)
    do i = imax - 3, imax - 2 ! DO i = imax-2,imax-2
        a(i) = 0.0_wp
        b(i) = -c0136
        c(i) = -c1418 + lambda*c0103*dx(i - 1)
        d(i) = 0.0_wp + lambda*dx(i)
        e(i) = c1418 + lambda*c0103*dx(i + 1)
        f(i) = c0136
        g(i) = 0.0_wp
    end do
! fifth-order biased
    a(imax - 1) = 0.0_wp
    b(imax - 1) = -1.0_wp/18.0_wp
    c(imax - 1) = -1.0_wp + lambda*1.0_wp/2.0_wp*dx(imax - 2)
    d(imax - 1) = 1.0_wp/2.0_wp + lambda*dx(imax - 1)
    e(imax - 1) = 5.0_wp/9.0_wp + lambda*1.0_wp/6.0_wp*dx(imax)
    f(imax - 1) = 0.0_wp
    g(imax - 1) = 0.0_wp
! third-order biased
    a(imax) = 0.0_wp
    b(imax) = -0.5_wp
    c(imax) = -2.0_wp + lambda*2.0_wp*dx(imax - 1)
    d(imax) = 5.0_wp/2.0_wp + lambda*dx(imax)
    e(imax) = 0.0_wp
    f(imax) = 0.0_wp
    g(imax) = 0.0_wp

! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
    if (ibc == 1) then
        d(1) = -5.0_wp/2.0_wp                                   ! array B22R
        d(2) = -5.0_wp/6.0_wp + lambda*2.0_wp/3.0_wp*dx(2)      ! B22R + lambda A22R fifth-order biased
        e(2) = 11.0_wp/12.0_wp + lambda*1.0_wp/2.0_wp*dx(3)
        f(imax) = 1.0_wp/6.0_wp                                   ! element A(imax-1,imax)

    else if (ibc == 2) then
        b(1) = 1.0_wp/6.0_wp                                   ! element A(2,1)
        c(imax - 1) = -11.0_wp/12.0_wp + lambda*1.0_wp/2.0_wp*dx(imax - 2) ! B11R + lambda A11R fifth-order biased
        d(imax - 1) = 5.0_wp/6.0_wp + lambda*2.0_wp/3.0_wp*dx(imax - 1)
        d(imax) = 5.0_wp/2.0_wp                                   ! array B11R

    end if

! -------------------------------------------------------------------
! Setting the RHS for the null space (with minus sign)
! -------------------------------------------------------------------
    h = 0.0_wp
    if (ibc == 1) then
        h(1) = 1.0_wp       ! normalization
        h(2) = c0136*5.0_wp ! fifth-order biased
        h(3) = c0136
    else if (ibc == 2) then
        h(imax - 2) = -c0136       ! fifth-order biased
        h(imax - 1) = -c0136*5.0_wp
        h(imax) = 1.0_wp       ! normalization
    end if

    return
end subroutine INT_C1N6M_LHS_E

!########################################################################
!Left-hand side case \lambda=0; pentadiagonal matrix of the linear system
!########################################################################
subroutine INT_C1N6M_LHS(imax, ibc, a, b, c, d, e, f, g)
    use TLab_Constants, only: wp, wi

    use FDM_PROCS, only: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6

    implicit none

    integer(wi), intent(IN) :: imax, ibc
    real(wp), dimension(imax), intent(OUT) :: a, b, c, d, e, f, g

! -------------------------------------------------------------------
    real(wp) :: c0136, c1418
    integer(wi) :: i

! ###################################################################
    c0136 = 1.0_wp/36.0_wp
    c1418 = 14.0_wp/18.0_wp

! -------------------------------------------------------------------
! Define diagonals of heptadiagonal system
! -------------------------------------------------------------------
! third-order biased
    a(1) = 0.0_wp
    b(1) = 0.0_wp
    c(1) = 0.0_wp
    d(1) = -5.0_wp/2.0_wp
    e(1) = 2.0_wp
    f(1) = 0.5_wp
    g(1) = 0.0_wp
! fifth-order biased
    a(2) = 0.0_wp
    b(2) = 0.0_wp
    c(2) = -5.0_wp/9.0_wp
    d(2) = -1.0_wp/2.0_wp
    e(2) = 1.0_wp
    f(2) = 1.0_wp/18.0_wp
    g(2) = 0.0_wp
! sixth-order centered (alpha=1/3)
    do i = 3, 4 ! DO i = 3,3
        a(i) = 0.0_wp
        b(i) = -c0136
        c(i) = -c1418
        d(i) = 0.0_wp
        e(i) = c1418
        f(i) = c0136
        g(i) = 0.0_wp
    end do
! sixth-order modified centered
    do i = 5, imax - 4 ! DO i = 4,imax-3
        a(i) = -C1N6M_CD6
        b(i) = -C1N6M_BD4
        c(i) = -C1N6M_AD2
        d(i) = 0.0_wp
        e(i) = C1N6M_AD2
        f(i) = C1N6M_BD4
        g(i) = C1N6M_CD6
    end do
! sixth-order centered (alpha=1/3)
    do i = imax - 3, imax - 2 ! DO i = imax-2,imax-2
        a(i) = 0.0_wp
        b(i) = -c0136
        c(i) = -c1418
        d(i) = 0.0_wp
        e(i) = c1418
        f(i) = c0136
        g(i) = 0.0_wp
    end do
! fifth-order biased
    a(imax - 1) = 0.0_wp
    b(imax - 1) = -1.0_wp/18.0_wp
    c(imax - 1) = -1.0_wp
    d(imax - 1) = 1.0_wp/2.0_wp
    e(imax - 1) = 5.0_wp/9.0_wp
    f(imax - 1) = 0.0_wp
    g(imax - 1) = 0.0_wp
! third-order biased
    a(imax) = 0.0_wp
    b(imax) = -0.5_wp
    c(imax) = -2.0_wp
    d(imax) = 5.0_wp/2.0_wp
    e(imax) = 0.0_wp
    f(imax) = 0.0_wp
    g(imax) = 0.0_wp

! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
    if (ibc == 1) then ! array B22R
        d(2) = -5.0_wp/6.0_wp  ! fifth-order biased
        e(2) = 11.0_wp/12.0_wp
        f(imax) = 1.0_wp/6.0_wp

    else if (ibc == 2) then ! array B11R
        c(imax - 1) = -11.0_wp/12.0_wp ! fifth-order biased
        d(imax - 1) = 5.0_wp/6.0_wp
        b(1) = 1.0_wp/6.0_wp

    end if

    return
end subroutine INT_C1N6M_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
subroutine INT_C1N6M_RHS(imax, jkmax, ibc, dx, h, l)
    use TLab_Constants, only: wp, wi

    use FDM_PROCS, only: C1N6M_ALPHA, C1N6M_BETA

    implicit none

    integer(wi), intent(IN) :: imax, jkmax, ibc
    real(wp), dimension(imax), intent(IN) :: dx
    real(wp), dimension(jkmax, imax), intent(IN) :: h
    real(wp), dimension(jkmax, imax), intent(OUT) :: l

! -------------------------------------------------------------------
    integer(wi) :: i
    real(wp) :: c0102, c0103, c0104, c0106, c0203

! ###################################################################
    c0102 = 1.0_wp/2.0_wp
    c0103 = 1.0_wp/3.0_wp
    c0104 = 1.0_wp/4.0_wp
    c0106 = 1.0_wp/6.0_wp
    c0203 = 2.0_wp/3.0_wp

! -------------------------------------------------------------------
! Boundary conditions
! -------------------------------------------------------------------
    if (ibc == 1) then ! array A22R
! BCs, see notes
        l(:, 1) = -2.0_wp*h(:, 2)*dx(2)                      ! contribution to u'_1
        l(:, 2) = c0102*h(:, 3)*dx(3) + c0203*h(:, 2)*dx(2) ! fifth-order biased
! fifth-order biased
        l(:, imax - 1) = c0102*h(:, imax - 2)*dx(imax - 2) + h(:, imax - 1)*dx(imax - 1) + c0106*h(:, imax)*dx(imax)
! third-order biased
        l(:, imax) = h(:, imax)*dx(imax) + 2.0_wp*h(:, imax - 1)*dx(imax - 1)

! -------------------------------------------------------------------
    else if (ibc == 2) then ! array A11R
! third-order biased
        l(:, 1) = h(:, 1)*dx(1) + 2.0_wp*h(:, 2)*dx(2)
! fifth-order
        l(:, 2) = c0106*h(:, 1)*dx(1) + h(:, 2)*dx(2) + c0102*h(:, 3)*dx(3)          ! fifth-order biased
! BCs, see notes
        l(:, imax - 1) = c0102*h(:, imax - 2)*dx(imax - 2) + c0203*h(:, imax - 1)*dx(imax - 1) ! fifth-order biased
        l(:, imax) = -2.0_wp*h(:, imax - 1)*dx(imax - 1)                                ! contribution to u'_N
    end if

! -------------------------------------------------------------------
! Interior points
! -------------------------------------------------------------------
! 6th-order centered with alpha=(1/3)
    do i = 3, 4 ! DO i = 3,3
        l(:, i) = h(:, i)*dx(i) + c0103*(h(:, i - 1)*dx(i - 1) + h(:, i + 1)*dx(i + 1))
    end do
    !
    do i = imax - 3, imax - 2 ! DO i = imax-2,imax-2
        l(:, i) = h(:, i)*dx(i) + c0103*(h(:, i - 1)*dx(i - 1) + h(:, i + 1)*dx(i + 1))
    end do

! sixth-order modified centered
    do i = 5, imax - 4 ! DO i = 4,imax-3
        l(:, i) = h(:, i)*dx(i) + &
                  C1N6M_ALPHA*(h(:, i - 1)*dx(i - 1) + h(:, i + 1)*dx(i + 1)) + &
                  C1N6M_BETA*(h(:, i - 2)*dx(i - 2) + h(:, i + 2)*dx(i + 2))
    end do

    return
end subroutine INT_C1N6M_RHS
