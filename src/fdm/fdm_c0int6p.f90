!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/22 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the finite difference interpolation with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, periodic.
!# Interior points according to Eq.C.1.4 (\alpha=3/10, \beta=0, c=0).
!# System multiplied by 4/3 to eliminate one multiplication in the RHS.
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
! #######################################################################
subroutine FDM_C0INT6P_LHS(imax, a, b, c)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax
    real(wp), dimension(imax), intent(OUT) :: a, b, c

! -------------------------------------------------------------------
    integer(wi) i

! ###################################################################
    do i = 1, imax
        a(i) = 2.0_wp/5.0_wp ! 3/10
        b(i) = 4.0_wp/3.0_wp ! 1
        c(i) = 2.0_wp/5.0_wp ! 3/10
    end do

    return
end subroutine FDM_C0INT6P_LHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from velocity to pressure grid
! #######################################################################
subroutine FDM_C0INTVP6P_RHS(imax, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, jkmax
    real(wp), dimension(jkmax, imax), intent(IN) :: u
    real(wp), dimension(jkmax, imax), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    integer(wi) :: im1, ip1, ip2, imm1
    real(wp) :: c0115

! #######################################################################

    c0115 = 1.0_wp/15.0_wp

    imm1 = imax - 1
    do i = 1, imax
        im1 = i - 1; im1 = im1 + imm1; im1 = MOD(im1, imax) + 1
        ip1 = i + 1; ip1 = ip1 + imm1; ip1 = MOD(ip1, imax) + 1
        ip2 = i + 2; ip2 = ip2 + imm1; ip2 = MOD(ip2, imax) + 1

        do jk = 1, jkmax
            d(jk, i) = u(jk, ip1) + u(jk, i) + c0115*(u(jk, ip2) + u(jk, im1))
        end do

    end do

    return
end subroutine FDM_C0INTVP6P_RHS

! #######################################################################
! Right-hand side; forcing term
! ==> interpolation from pressure to velocity grid
! #######################################################################
subroutine FDM_C0INTPV6P_RHS(imax, jkmax, u, d)
    use TLab_Constants, only: wp, wi

    implicit none

    integer(wi), intent(IN) :: imax, jkmax
    real(wp), dimension(jkmax, imax), intent(IN) :: u
    real(wp), dimension(jkmax, imax), intent(OUT) :: d

! -------------------------------------------------------------------
    integer(wi) :: i, jk
    integer(wi) :: im1, ip1, im2, imm1
    real(wp) :: c0115

! #######################################################################

    c0115 = 1.0_wp/15.0_wp

    imm1 = imax - 1
    do i = 1, imax
        im1 = i - 1; im1 = im1 + imm1; im1 = MOD(im1, imax) + 1
        im2 = i - 2; im2 = im2 + imm1; im2 = MOD(im2, imax) + 1
        ip1 = i + 1; ip1 = ip1 + imm1; ip1 = MOD(ip1, imax) + 1

        do jk = 1, jkmax
            d(jk, i) = u(jk, i) + u(jk, im1) + c0115*(u(jk, ip1) + u(jk, im2))
        end do

    end do

    return
end subroutine FDM_C0INTPV6P_RHS
