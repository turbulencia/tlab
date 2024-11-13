!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2022/01/14 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the finite difference interpolation with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, non-periodic.
!# Interior points according to Eq.C.1.4 (\alpha=3/10, \beta=0, c=0).
!#
!# Different boundary closures can be found in Albin 2010
!# (https://doi.org/10.1002/fld.2520) table 6 (typos in Si2_5i and Si2_6i).
!#
!########################################################################
!# ARGUMENTS
!#
!# u    In    function to be interpolated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! ==> interpolation from velocity to pressure grid
! #######################################################################
subroutine FDM_C0INTVP6_LHS(imaxp, a, b, c)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: imaxp ! pressure grid (imaxp==imax-1)
    real(wp), dimension(imaxp), intent(OUT) :: a, b, c

! -------------------------------------------------------------------
    integer(wi) :: i

! ###################################################################
! forth-order biased (implicit)
    a(1) = 0.0_wp
    b(1) = 1.0_wp
    c(1) = 1.0_wp
    !
    a(imaxp) = 1.0_wp
    b(imaxp) = 1.0_wp
    c(imaxp) = 0.0_wp

! sixth-order (implicit)
    do i = 2, imaxp - 1
        a(i) = 3.0_wp/10.0_wp
        b(i) = 1.0_wp
        c(i) = 3.0_wp/10.0_wp
    end do

! ! forth-order (implicit) - for adjacent boundary nodes
!   a(2)      = 1.0_wp / 6.0_wp
!   b(2)      = 1.0_wp
!   c(2)      = 1.0_wp / 6.0_wp
!   !
!   a(imaxp-1) = 1.0_wp / 6.0_wp
!   b(imaxp-1) = 1.0_wp
!   c(imaxp-1) = 1.0_wp / 6.0_wp

    return
end subroutine FDM_C0INTVP6_LHS

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! ==> interpolation from pressure to velocity grid
! #######################################################################
subroutine FDM_C0INTPV6_LHS(imax, a, b, c)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax ! velocity grid
    real(wp), dimension(imax), intent(OUT) :: a, b, c

! -------------------------------------------------------------------
    integer(wi) :: i

! ###################################################################
! third-order biased (implicit)
    a(1) = 0.0_wp
    b(1) = 1.0_wp
    c(1) = 3.0_wp
    !
    a(imax) = 3.0_wp
    b(imax) = 1.0_wp
    c(imax) = 0.0_wp

! forth-order (implicit)
    a(2) = 1.0_wp/6.0_wp
    b(2) = 1.0_wp
    c(2) = 1.0_wp/6.0_wp
    !
    a(imax - 1) = 1.0_wp/6.0_wp
    b(imax - 1) = 1.0_wp
    c(imax - 1) = 1.0_wp/6.0_wp

! sixth-order (implicit)
    do i = 3, imax - 2
        a(i) = 3.0_wp/10.0_wp
        b(i) = 1.0_wp
        c(i) = 3.0_wp/10.0_wp
    end do

    return
end subroutine FDM_C0INTPV6_LHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from velocity to pressure grid
! #######################################################################
subroutine FDM_C0INTVP6_RHS(imax, imaxp, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, imaxp, jkmax
    real(wp), dimension(jkmax, imax), intent(IN) :: u
    real(wp), dimension(jkmax, imaxp), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    real(wp) :: c0302, c0104, c0304, c0120 !, c0203

! #######################################################################
    c0302 = 3.0_wp/2.0_wp
    c0104 = 1.0_wp/4.0_wp
    c0304 = 3.0_wp/4.0_wp
    c0120 = 1.0_wp/20.0_wp
    ! c0203 = 2.0_wp / 3.0_wp

    do jk = 1, jkmax
        ! forth-order biased (implicit)
        d(jk, 1) = c0302*u(jk, 2) + c0104*(u(jk, 1) + u(jk, 3))
        d(jk, imaxp) = c0302*u(jk, imax - 1) + c0104*(u(jk, imax) + u(jk, imax - 2))
        ! forth-order (implicit)
        ! d(jk,2)       = c0203*(u(jk,3)      + u(jk,2)     )
        ! d(jk,imaxp-1) = c0203*(u(jk,imax-2) + u(jk,imax-1))
    end do

! sixth-order (implicit)
    do i = 2, imaxp - 1
        ! DO i = 3,imaxp-2 ! for forth order (implict)
        do jk = 1, jkmax
            d(jk, i) = c0304*(u(jk, i + 1) + u(jk, i)) + c0120*(u(jk, i + 2) + u(jk, i - 1))
        end do
    end do

    return
end subroutine FDM_C0INTVP6_RHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from pressure to velocity grid
! #######################################################################
subroutine FDM_C0INTPV6_RHS(imax, imaxp, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, imaxp, jkmax
    real(wp), dimension(jkmax, imaxp), intent(IN) :: u
    real(wp), dimension(jkmax, imax), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    real(wp) :: c0304, c0120, c0203

! #######################################################################
    c0304 = 3.0_wp/4.0_wp
    c0120 = 1.0_wp/20.0_wp
    c0203 = 2.0_wp/3.0_wp

    do jk = 1, jkmax
        ! forth-order biased (implicit)
        d(jk, 1) = 3.0_wp*u(jk, 1) + u(jk, 2)
        d(jk, imax) = 3.0_wp*u(jk, imaxp) + u(jk, imaxp - 1)
        ! forth-order (implicit)
        d(jk, 2) = c0203*(u(jk, 1) + u(jk, 2))
        d(jk, imax - 1) = c0203*(u(jk, imaxp) + u(jk, imaxp - 1))
    end do

! sixth-order (implicit)
    do i = 3, imax - 2
        do jk = 1, jkmax
            d(jk, i) = c0304*(u(jk, i) + u(jk, i - 1)) + c0120*(u(jk, i + 1) + u(jk, i - 2))
        end do
    end do

    return
end subroutine FDM_C0INTPV6_RHS
