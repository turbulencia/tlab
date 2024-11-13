!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/01/14 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the interpolatory finite difference  first derivative
!# with 6th-order tridiagonal compact scheme by JCP Lele 1992, non-periodic.
!# Interior points according to Eq. B.1.1 (\alpha=9/62, \beta=0, c=0).
!# System for this scheme is multiplied by 62/63 to eliminate one
!# multiplication in the RHS.
!#
!# Different boundary closures can be found in Albin 2010
!# (https://doi.org/10.1002/fld.2520) table 6 (typo in Sd2_3i).
!#
!########################################################################
!# ARGUMENTS
!#
!# u    In    function to be diferentiated
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! ==> interpolation from velocity to pressure grid
! #######################################################################
subroutine FDM_C1INTVP6_LHS(imaxp, dx, a, b, c)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imaxp ! pressure grid
    real(wp), dimension(imaxp), intent(IN) :: dx    ! pressure grid
    real(wp), dimension(imaxp), intent(OUT) :: a, b, c

! -------------------------------------------------------------------
    integer(wi) :: i

! ###################################################################
! ! third-order biased (explicit)
!   a(1)     = 0.0_wp
!   b(1)     = 1.0_wp
!   c(1)     = 0.0_wp
!   !
!   a(imaxp) = 0.0_wp
!   b(imaxp) = 1.0_wp
!   c(imaxp) = 0.0_wp

! third-order biased (implicit)
    a(1) = 0.0_wp
    b(1) = 1.0_wp
    c(1) = -1.0_wp
    !
    a(imaxp) = -1.0_wp
    b(imaxp) = 1.0_wp
    c(imaxp) = 0.0_wp

! sixth-order (implicit)
    do i = 2, imaxp - 1
        a(i) = 9.0_wp/63.0_wp ! 9/62
        b(i) = 62.0_wp/63.0_wp ! 1
        c(i) = 9.0_wp/63.0_wp ! 9/62
    end do

! ! forth-order (implicit) - for adjacent boundary nodes
!   a(2)       = 1.0_wp / 22.0_wp
!   b(2)       = 1.0_wp
!   c(2)       = 1.0_wp / 22.0_wp
!   !
!   a(imaxp-1) = 1.0_wp / 22.0_wp
!   b(imaxp-1) = 1.0_wp
!   c(imaxp-1) = 1.0_wp / 22.0_wp

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
    c(imaxp) = c(imaxp)*dx(1)
    b(1) = b(1)*dx(1)
    a(2) = a(2)*dx(1)

    do i = 2, imaxp - 1
        c(i - 1) = c(i - 1)*dx(i)
        b(i) = b(i)*dx(i)
        a(i + 1) = a(i + 1)*dx(i)
    end do

    c(imaxp - 1) = c(imaxp - 1)*dx(imaxp)
    b(imaxp) = b(imaxp)*dx(imaxp)
    a(1) = a(1)*dx(imaxp)

    return
end subroutine FDM_C1INTVP6_LHS

! #######################################################################
! Left-hand side; tridiagonal matrix of the linear system
! ==> interpolation from pressure to velocity grid
! #######################################################################
subroutine FDM_C1INTPV6_LHS(imax, dx, a, b, c)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax  ! velocity grid
    real(wp), dimension(imax), intent(IN) :: dx    ! veloctiy grid
    real(wp), dimension(imax), intent(OUT) :: a, b, c

! -------------------------------------------------------------------
    integer(wi) :: i

! ###################################################################
! ! third-order biased (explicit)
!   a(1)    = 0.0_wp
!   b(1)    = 1.0_wp
!   c(1)    = 0.0_wp
!   !
!   a(imax) = 0.0_wp
!   b(imax) = 1.0_wp
!   c(imax) = 0.0_wp

! third-order biased (implicit)
    a(1) = 0.0_wp
    b(1) = 1.0_wp
    c(1) = 23.0_wp
    !
    a(imax) = 23.0_wp
    b(imax) = 1.0_wp
    c(imax) = 0.0_wp

    ! forth-order (implicit) - for adjacent boundary nodes
    a(2) = 1.0_wp/22.0_wp
    b(2) = 1.0_wp
    c(2) = 1.0_wp/22.0_wp
    !
    a(imax - 1) = 1.0_wp/22.0_wp
    b(imax - 1) = 1.0_wp
    c(imax - 1) = 1.0_wp/22.0_wp

! sixth-order (implicit)
    do i = 3, imax - 2
        a(i) = 9.0_wp/63.0_wp ! 9/62
        b(i) = 62.0_wp/63.0_wp ! 1
        c(i) = 9.0_wp/63.0_wp ! 9/62
    end do

! -------------------------------------------------------------------
! Jacobian Multiplication
! -------------------------------------------------------------------
    c(imax) = c(imax)*dx(1)
    b(1) = b(1)*dx(1)
    a(2) = a(2)*dx(1)

    do i = 2, imax - 1
        c(i - 1) = c(i - 1)*dx(i)
        b(i) = b(i)*dx(i)
        a(i + 1) = a(i + 1)*dx(i)
    end do

    c(imax - 1) = c(imax - 1)*dx(imax)
    b(imax) = b(imax)*dx(imax)
    a(1) = a(1)*dx(imax)

    return
end subroutine FDM_C1INTPV6_LHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from velocity to pressure grid
! #######################################################################
subroutine FDM_C1INTVP6_RHS(imax, imaxp, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, imaxp, jkmax
    real(wp), dimension(jkmax, imax), intent(IN) :: u
    real(wp), dimension(jkmax, imaxp), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    real(wp) :: c17189, c0120
    ! real(wp)                                        :: c0108, c0708, c0124, c2324, c1211

! #######################################################################
    c17189 = 17.0_wp/189.0_wp
    c0120 = 1.0_wp/20.0_wp

    ! c0108  = 1.0_wp  / 8.0_wp
    ! c0708  = C_7_R  / 8.0_wp
    ! c0124  = 1.0_wp  / C_24_R
    ! c2324  = 23.0_wp / C_24_R
    ! c1211  = 12.0_wp / 11.0_wp

    do jk = 1, jkmax
        ! ! third-order biased (explicit)
        ! d(jk,1)       = - c2324*u(jk,1)    + c0708*u(jk,2)      + c0108*u(jk,3)      - c0124*u(jk,4)
        ! d(jk,imaxp)   =   c2324*u(jk,imax) - c0708*u(jk,imax-1) - c0108*u(jk,imax-2) + c0124*u(jk,imax-3)
        ! forth-order biased (implicit)
        d(jk, 1) = -u(jk, 1) + 2.0_wp*u(jk, 2) - u(jk, 3)
        d(jk, imaxp) = u(jk, imax) - 2.0_wp*u(jk, imax - 1) + u(jk, imax - 2)
        ! ! forth-order (implicit)
        ! d(jk,2)       =   c1211 * (u(jk,3)       - u(jk,2))
        ! d(jk,imaxp-1) = - c1211 * (u(jk,imax-2) - u(jk,imax-1))
    end do

    ! sixth-order (implicit)
    do i = 2, imaxp - 1
        ! DO i = 3,imaxp-2 ! for forth order (implict)
        do jk = 1, jkmax
            d(jk, i) = (u(jk, i + 1) - u(jk, i)) + c17189*(u(jk, i + 2) - u(jk, i - 1))
        end do
    end do

    return
end subroutine FDM_C1INTVP6_RHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from pressure to velocity grid
! #######################################################################
subroutine FDM_C1INTPV6_RHS(imax, imaxp, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, imaxp, jkmax
    real(wp), dimension(jkmax, imaxp), intent(IN) :: u
    real(wp), dimension(jkmax, imax), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    real(wp) :: c17189, c0120, c1211
    ! real(wp)                                        :: c0108, c0708, c0124, c2324

! #######################################################################
    c17189 = 17.0_wp/189.0_wp
    c0120 = 1.0_wp/20.0_wp
    c1211 = 12.0_wp/11.0_wp

    ! c2324  = 23.0_wp / C_24_R
    ! c0708  = C_7_R  / 8.0_wp
    ! c0108  = 1.0_wp  / 8.0_wp
    ! c0124  = 1.0_wp  / C_24_R

    do jk = 1, jkmax
        ! ! third-order biased (explicit)
        ! d(jk,1)      = - c2324*u(jk,1)     + c0708*u(jk,2)       + c0108*u(jk,3)       - c0124*u(jk,4)
        ! d(jk,imax)   =   c2324*u(jk,imaxp) - c0708*u(jk,imaxp-1) - c0108*u(jk,imaxp-2) + c0124*u(jk,imaxp-3)
        ! forth-order biased (implicit)
        d(jk, 1) = -25.0_wp*u(jk, 1) + 26.0_wp*u(jk, 2) - u(jk, 3)
        d(jk, imax) = 25.0_wp*u(jk, imaxp) - 26.0_wp*u(jk, imaxp - 1) + u(jk, imaxp - 2)
        ! forth-order (implicit)
        d(jk, 2) = c1211*(u(jk, 2) - u(jk, 1))
        d(jk, imax - 1) = -c1211*(u(jk, imaxp - 1) - u(jk, imaxp))
    end do

    ! sixth-order (implicit)
    do i = 3, imax - 2
        ! DO i = 3,imaxp-2 ! for forth order (implict)
        do jk = 1, jkmax
            d(jk, i) = (u(jk, i) - u(jk, i - 1)) + c17189*(u(jk, i + 1) - u(jk, i - 2))
        end do
    end do

    return
end subroutine FDM_C1INTPV6_RHS
