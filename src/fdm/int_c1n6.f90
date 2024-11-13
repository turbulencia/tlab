!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2010/10/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the first derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic,
!# in order to solve the IVP
!#
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# The system of N-1 eqns:
!#
!#                    (B + \lambda A)u = Af = g
!#
!# is established in this routine, giving diagonals a-e and g (see notes).
!# Interior points 6th-order according to Eq. 2.1.7.
!# The second point from Eq. 2.1.6 forth-order (b=0).
!# The first point from third-order biased Eq. 4.1.3 (d=0).
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
!# f           Out   forcing term for the exponential
!#
!########################################################################

!########################################################################
!Left-hand side; pentadiagonal matrix of the linear system and f
!########################################################################
subroutine INT_C1N6_LHS_E(imax, ibc, dx, lambda, a, b, c, d, e, f)
    use TLab_Constants, only: wp, wi

    implicit none

    real(wp) lambda
    integer(wi), intent(IN) :: imax, ibc
    real(wp), dimension(imax), intent(IN) :: dx
    real(wp), dimension(imax), intent(OUT) :: a, b, c, d, e, f

! -------------------------------------------------------------------
    integer(wi) i
    real(wp) c0136, c1418, c0103, c0104

! ###################################################################
    c0136 = 1.0_wp/36.0_wp
    c1418 = 14.0_wp/18.0_wp
    c0103 = 1.0_wp/3.0_wp
    c0104 = 1.0_wp/4.0_wp

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
    a(1) = 0.0_wp ! padding
    b(1) = 0.0_wp ! padding
    c(1) = -5.0_wp/2.0_wp + lambda*dx(1)
    d(1) = 2.0_wp + lambda*2.0_wp*dx(2)
    e(1) = 0.5_wp
! fourth-order centered
!   a(2       ) = 0.0_wp ! padding
!   b(2       ) =-3.0_wp/4.0_wp + lambda*c0104 *dx(1)
!   c(2       ) = 0.0_wp       + lambda       *dx(2)
!   d(2       ) = 3.0_wp/4.0_wp + lambda*c0104 *dx(3)
!   e(2       ) = 0.0_wp
! fifth-order biased
    a(2) = 0.0_wp ! padding
    b(2) = -5.0_wp/9.0_wp + lambda*1.0_wp/6.0_wp*dx(1)
    c(2) = -1.0_wp/2.0_wp + lambda*dx(2)
    d(2) = 1.0_wp + lambda*1.0_wp/2.0_wp*dx(3)
    e(2) = 1.0_wp/18.0_wp
! sixth-order centered
    do i = 3, imax - 2
        a(i) = -c0136
        b(i) = -c1418 + lambda*c0103*dx(i - 1)
        c(i) = 0.0_wp + lambda*dx(i)
        d(i) = c1418 + lambda*c0103*dx(i + 1)
        e(i) = c0136
    end do
! fourth-order centered
!   a(imax-1  ) = 0.0_wp
!   b(imax-1  ) =-3.0_wp/4.0_wp + lambda*c0104 *dx(imax-2)
!   c(imax-1  ) = 0.0_wp       + lambda       *dx(imax-1)
!   d(imax-1  ) = 3.0_wp/4.0_wp + lambda*c0104 *dx(imax  )
!   e(imax-1  ) = 0.0_wp ! padding
! fifth-order biased
    a(imax - 1) = -1.0_wp/18.0_wp
    b(imax - 1) = -1.0_wp + lambda*1.0_wp/2.0_wp*dx(imax - 2)
    c(imax - 1) = 1.0_wp/2.0_wp + lambda*dx(imax - 1)
    d(imax - 1) = 5.0_wp/9.0_wp + lambda*1.0_wp/6.0_wp*dx(imax)
    e(imax - 1) = 0.0_wp ! padding
! third-order biased
    a(imax) = -0.5_wp
    b(imax) = -2.0_wp + lambda*2.0_wp*dx(imax - 1)
    c(imax) = 5.0_wp/2.0_wp + lambda*dx(imax)
    d(imax) = 0.0_wp ! padding
    e(imax) = 0.0_wp ! padding

! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
    if (ibc == 1) then
        c(1) = -5.0_wp/2.0_wp                            ! array B22R
!     c(2     ) =-1.0_wp/2.0_wp + lambda*1.0_wp/2.0_wp*dx(2) ! array B22R + lambda A22R
!     d(2     ) = 5.0_wp/8.0_wp + lambda*1.0_wp/4.0_wp*dx(3)
!     e(imax  ) = C_025_R                                 ! element A(imax-1,imax)
        c(2) = -5.0_wp/6.0_wp + lambda*2.0_wp/3.0_wp*dx(2) ! B22R + lambda A22R fifth-order biased
        d(2) = 11.0_wp/12.0_wp + lambda*1.0_wp/2.0_wp*dx(3)
        e(imax) = 1.0_wp/6.0_wp                              ! element A(imax-1,imax)

    else if (ibc == 2) then
        a(1) = 1.0_wp/6.0_wp                                   ! element A(2,1)
        b(imax - 1) = -11.0_wp/12.0_wp + lambda*1.0_wp/2.0_wp*dx(imax - 2) ! B11R + lambda A11R fifth-order biased
        c(imax - 1) = 5.0_wp/6.0_wp + lambda*2.0_wp/3.0_wp*dx(imax - 1)
!     a(1)      = C_025_R                                     ! element A(2,1) 4th
!     b(imax-1) =-5.0_wp/8.0_wp + lambda*1.0_wp/4.0_wp*dx(imax-2) ! array B11R + lambda A11R 4th
!     c(imax-1) = 1.0_wp/2.0_wp + lambda*1.0_wp/2.0_wp*dx(imax-1)
        c(imax) = 5.0_wp/2.0_wp                                 ! array B11R

    end if

! -------------------------------------------------------------------
! Setting the RHS for the null space (with minus sign)
! -------------------------------------------------------------------
    f = 0.0_wp
    if (ibc == 1) then
        f(1) = 1.0_wp       ! normalization
!     f(2     ) = 1.0_wp/8.0_wp ! b21R
!     f(3     ) = c0136
        f(2) = c0136*5.0_wp ! fifth-order biased
        f(3) = c0136
    else if (ibc == 2) then
!     f(imax-2) =-c0136       ! b1NR
!     f(imax-1) =-1.0_wp/8.0_wp
        f(imax - 2) = -c0136       ! fifth-order biased
        f(imax - 1) = -c0136*5.0_wp
        f(imax) = 1.0_wp       ! normalization
    end if

    return
end subroutine INT_C1N6_LHS_E

!########################################################################
!Left-hand side case \lambda=0; pentadiagonal matrix of the linear system
!########################################################################
subroutine INT_C1N6_LHS(imax, ibc, a, b, c, d, e)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, ibc
    real(wp), dimension(imax), intent(OUT) :: a, b, c, d, e

! -------------------------------------------------------------------
    real(wp) c0136, c1418

! ###################################################################
    c0136 = 1.0_wp/36.0_wp
    c1418 = 14.0_wp/18.0_wp

! -------------------------------------------------------------------
! Define diagonals of pentadiagonal system
! -------------------------------------------------------------------
! third-order biased
    a(1) = 0.0_wp ! padding
    b(1) = 0.0_wp ! padding
    c(1) = -5.0_wp/2.0_wp
    d(1) = 2.0_wp
    e(1) = 0.5_wp
! fourth-order centered
!    a(2       ) = 0.0_wp ! padding
!    b(2       ) =-3.0_wp/4.0_wp
!    c(2       ) = 0.0_wp
!    d(2       ) = 3.0_wp/4.0_wp
!    e(2       ) = 0.0_wp
! fifth-order biased
    a(2) = 0.0_wp ! padding
    b(2) = -5.0_wp/9.0_wp
    c(2) = -1.0_wp/2.0_wp
    d(2) = 1.0_wp
    e(2) = 1.0_wp/18.0_wp
! sixth-order centered
    a(3:imax - 2) = -c0136
    b(3:imax - 2) = -c1418
    c(3:imax - 2) = 0.0_wp
    d(3:imax - 2) = c1418
    e(3:imax - 2) = c0136
! fourth-order centered
!   a(  imax-1) = 0.0_wp
!   b(  imax-1) =-3.0_wp/4.0_wp
!   c(  imax-1) = 0.0_wp
!   d(  imax-1) = 3.0_wp/4.0_wp
!   e(  imax-1) = 0.0_wp ! padding
! fifth-order biased
    a(imax - 1) = -1.0_wp/18.0_wp
    b(imax - 1) = -1.0_wp
    c(imax - 1) = 1.0_wp/2.0_wp
    d(imax - 1) = 5.0_wp/9.0_wp
    e(imax - 1) = 0.0_wp ! padding
! third-order biased
    a(imax) = -0.5_wp
    b(imax) = -2.0_wp
    c(imax) = 5.0_wp/2.0_wp
    d(imax) = 0.0_wp ! padding
    e(imax) = 0.0_wp ! padding

! -------------------------------------------------------------------
! Boundary conditions, see notes
! -------------------------------------------------------------------
    if (ibc == 1) then ! array B22R
!     c(2     ) =-1.0_wp/2.0_wp ! fourth-order biased
!     d(2     ) = 5.0_wp/8.0_wp
!     e(imax  ) = C_025_R     ! element A(imax-1,imax)
        c(2) = -5.0_wp/6.0_wp ! fifth-order biased
        d(2) = 11.0_wp/12.0_wp
        e(imax) = 1.0_wp/6.0_wp

    else if (ibc == 2) then ! array B11R
!     b(imax-1) =-5.0_wp/8.0_wp ! fourth-order biased
!     c(imax-1) = 1.0_wp/2.0_wp
!     a(1)      = C_025_R     ! element A(2,1)
        b(imax - 1) = -11.0_wp/12.0_wp ! fifth-order biased
        c(imax - 1) = 5.0_wp/6.0_wp
        a(1) = 1.0_wp/6.0_wp

    end if

    return
end subroutine INT_C1N6_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
subroutine INT_C1N6_RHS(imax, jkmax, ibc, dx, f, g)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, jkmax, ibc
    real(wp), dimension(imax), intent(IN) :: dx
    real(wp), dimension(jkmax, imax), intent(IN) :: f
    real(wp), dimension(jkmax, imax), intent(OUT) :: g

! -------------------------------------------------------------------
    integer(wi) i
    real(wp) c0102, c0103, c0104, c0106, c0203

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
        g(:, 1) = -2.0_wp*f(:, 2)*dx(2) ! contribution to u'_1
!     g(:,2) = c0104*f(:,3)*dx(3) + c0102*f(:,2)*dx(2)
        g(:, 2) = c0102*f(:, 3)*dx(3) + c0203*f(:, 2)*dx(2) ! fifth-order biased

! fourth-order centered
!     g(:,imax-1) = c0104*(f(:,imax-2)*dx(imax-2)+f(:,imax)*dx(imax)) + f(:,imax-1)*dx(imax-1)
! fifth-order biased
        g(:, imax - 1) = c0102*f(:, imax - 2)*dx(imax - 2) + f(:, imax - 1)*dx(imax - 1) + c0106*f(:, imax)*dx(imax)
! third-order biased
        g(:, imax) = f(:, imax)*dx(imax) + 2.0_wp*f(:, imax - 1)*dx(imax - 1)

! -------------------------------------------------------------------
    else if (ibc == 2) then ! array A11R
! third-order biased
        g(:, 1) = f(:, 1)*dx(1) + 2.0_wp*f(:, 2)*dx(2)
! fourth-order centered
!     g(:,2) = c0104*(f(:,1)*dx(1)+f(:,3)*dx(3)) + f(:,2)*dx(2)
! fifth-order
        g(:, 2) = c0106*f(:, 1)*dx(1) + f(:, 2)*dx(2) + c0102*f(:, 3)*dx(3) ! fifth-order biased

! BCs, see notes
!     g(:,imax-1) = c0104*f(:,imax-2)*dx(imax-2) + c0102*f(:,imax-1)*dx(imax-1)
        g(:, imax - 1) = c0102*f(:, imax - 2)*dx(imax - 2) + c0203*f(:, imax - 1)*dx(imax - 1)! fifth-order biased
        g(:, imax) = -2.0_wp*f(:, imax - 1)*dx(imax - 1) ! contribution to u'_N
    end if

! -------------------------------------------------------------------
! Interior points
! -------------------------------------------------------------------
! sixth-order
    do i = 3, imax - 2
        g(:, i) = f(:, i)*dx(i) + c0103*(f(:, i - 1)*dx(i - 1) + f(:, i + 1)*dx(i + 1))
    end do

    return
end subroutine INT_C1N6_RHS
