!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/23 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the interpolatory finite difference  first derivative
!# with 6th-order tridiagonal compact scheme by JCP Lele 1992, periodic.
!# Interior points according to Eq. B.1.1 (\alpha=9/62, \beta=0, c=0).
!# System multiplied by 62/63 to eliminate one multiplication in the RHS.
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
! #######################################################################
subroutine FDM_C1INT6P_LHS(imax, dx, a, b, c)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax
    real(wp), dimension(imax), intent(IN) :: dx
    real(wp), dimension(imax), intent(OUT) :: a, b, c

! -------------------------------------------------------------------
    integer(wi) i

! ###################################################################
    do i = 1, imax
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
end subroutine FDM_C1INT6P_LHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from velocity to pressure grid
! #######################################################################
subroutine FDM_C1INTVP6P_RHS(imax, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, jkmax
    real(wp), dimension(jkmax, imax), intent(IN) :: u
    real(wp), dimension(jkmax, imax), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    integer(wi) :: im1, ip1, ip2, imm1
    real(wp) :: c17189

! #######################################################################

    c17189 = 17.0_wp/189.0_wp

    imm1 = imax - 1
    do i = 1, imax
        im1 = i - 1; im1 = im1 + imm1; im1 = MOD(im1, imax) + 1
        ip1 = i + 1; ip1 = ip1 + imm1; ip1 = MOD(ip1, imax) + 1
        ip2 = i + 2; ip2 = ip2 + imm1; ip2 = MOD(ip2, imax) + 1

        do jk = 1, jkmax
            d(jk, i) = (u(jk, ip1) - u(jk, i)) + c17189*(u(jk, ip2) - u(jk, im1))
        end do

    end do

    return
end subroutine FDM_C1INTVP6P_RHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from pressure to velocity grid
! #######################################################################
subroutine FDM_C1INTPV6P_RHS(imax, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, jkmax
    real(wp), dimension(jkmax, imax), intent(IN) :: u
    real(wp), dimension(jkmax, imax), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    integer(wi) :: im1, ip1, im2, imm1
    real(wp) :: c17189

! #######################################################################

    c17189 = 17.0_wp/189.0_wp

    imm1 = imax - 1
    do i = 1, imax
        im1 = i - 1; im1 = im1 + imm1; im1 = MOD(im1, imax) + 1
        im2 = i - 2; im2 = im2 + imm1; im2 = MOD(im2, imax) + 1
        ip1 = i + 1; ip1 = ip1 + imm1; ip1 = MOD(ip1, imax) + 1

        do jk = 1, jkmax
            d(jk, i) = (u(jk, i) - u(jk, im1)) + c17189*(u(jk, ip1) - u(jk, im2))
        end do

    end do

    return
end subroutine FDM_C1INTPV6P_RHS
