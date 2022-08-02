#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# 4th-Order Compact Filter from Lele [J. Comp. Phys., V103, 1992]
!#
!# uf_i + alpha*(uf_i-1 + uf_i+1) = a*u_i + b*(u_i-1 + u_i+1) + c*(u_i-2 + u_i+2)
!#
!########################################################################
subroutine FLT_C4_LHS(imax, bcsimin, bcsimax, alpha, a,b,c)
    implicit none
    TINTEGER, intent(in) :: imax, bcsimin, bcsimax
    TREAL, intent(in) :: alpha
    TREAL, intent(out) :: a(imax), b(imax), c(imax)

    a(:) = alpha
    b(:) = C_1_R    
    c(:) = alpha

    if ( bcsimin == DNS_FILTER_BCS_ZERO ) then
        c(1) = C_0_R
    end if

    if ( bcsimax== DNS_FILTER_BCS_ZERO ) then
        a(imax) = C_0_R
    end if

    return
end subroutine FLT_C4_LHS

subroutine FLT_C4_RHS(imax, jkmax, periodic, bcsimin, bcsimax, cxi, u, uf)!, txi)

    implicit none

    logical periodic
    TINTEGER, intent(IN) :: imax, jkmax, bcsimin, bcsimax
    TREAL, dimension(jkmax, imax), intent(IN) :: u
    TREAL, dimension(jkmax, imax), intent(OUT) :: uf
    TREAL, dimension(imax, 5), intent(IN) :: cxi

! -----------------------------------------------------------------------
    TINTEGER jk, i

! #######################################################################
! Set up right hand side and DO forward/backward substitution
! #######################################################################
    if (periodic) then ! periodic
        do jk = 1, jkmax
            uf(jk, 1) = cxi(1, 1)*u(jk, imax - 1) + cxi(1, 2)*u(jk, imax) + cxi(1, 3)*u(jk, 1) + &
                        cxi(1, 4)*u(jk, 2) + cxi(1, 5)*u(jk, 3)

            uf(jk, 2) = cxi(2, 1)*u(jk, imax) + cxi(2, 2)*u(jk, 1) + cxi(2, 3)*u(jk, 2) + &
                        cxi(2, 4)*u(jk, 3) + cxi(2, 5)*u(jk, 4)

            uf(jk, imax) = cxi(imax, 1)*u(jk, imax - 2) + cxi(imax, 2)*u(jk, imax - 1) &
                           + cxi(imax, 3)*u(jk, imax) + cxi(imax, 4)*u(jk, 1) + cxi(imax, 5)*u(jk, 2)

            uf(jk, imax - 1) = cxi(imax - 1, 1)*u(jk, imax - 3) + cxi(imax - 1, 2)*u(jk, imax - 2) &
                               + cxi(imax - 1, 3)*u(jk, imax - 1) + cxi(imax - 1, 4)*u(jk, imax) + cxi(imax - 1, 5)*u(jk, 1)
        end do

    else ! biased
        do jk = 1, jkmax
            uf(jk, 1) = cxi(1, 1)*u(jk, 1) &
                        + cxi(1, 2)*u(jk, 2) &
                        + cxi(1, 3)*u(jk, 3) &
                        + cxi(1, 4)*u(jk, 4) &
                        + cxi(1, 5)*u(jk, 5)
            uf(jk, 2) = cxi(2, 1)*u(jk, 1) &
                        + cxi(2, 2)*u(jk, 2) &
                        + cxi(2, 3)*u(jk, 3) &
                        + cxi(2, 4)*u(jk, 4) &
                        + cxi(2, 5)*u(jk, 5)

            uf(jk, imax - 1) = cxi(imax - 1, 5)*u(jk, imax) &
                               + cxi(imax - 1, 4)*u(jk, imax - 1) &
                               + cxi(imax - 1, 3)*u(jk, imax - 2) &
                               + cxi(imax - 1, 2)*u(jk, imax - 3) &
                               + cxi(imax - 1, 1)*u(jk, imax - 4)
            uf(jk, imax) = cxi(imax, 5)*u(jk, imax) &
                           + cxi(imax, 4)*u(jk, imax - 1) &
                           + cxi(imax, 3)*u(jk, imax - 2) &
                           + cxi(imax, 2)*u(jk, imax - 3) &
                           + cxi(imax, 1)*u(jk, imax - 4)
        end do

        if (bcsimin == DNS_FILTER_BCS_ZERO) then ! No filter at i=1
            do jk = 1, jkmax
                uf(jk, 1) = u(jk, 1)
            end do
        end if
        if (bcsimax == DNS_FILTER_BCS_ZERO) then ! No filter at i=1
            do jk = 1, jkmax
                uf(jk, imax) = u(jk, imax)
            end do
        end if

    end if

    do i = 3, imax - 2
        do jk = 1, jkmax
            uf(jk, i) = cxi(i, 1)*u(jk, i - 2) + cxi(i, 2)*u(jk, i - 1) + cxi(i, 3)*u(jk, i) + &
                        cxi(i, 4)*u(jk, i + 1) + cxi(i, 5)*u(jk, i + 2)
        end do
    end do

    return
end subroutine FLT_C4_RHS

subroutine FLT_C4P_CUTOFF_LHS(imax, a, b, c, d, e)
    implicit none

    TINTEGER, intent(in) :: imax
    TREAL, intent(out) :: a(imax), b(imax), c(imax), d(imax), e(imax)

    TREAL, parameter :: C4_ALPHA = 0.6522474d0
    TREAL, parameter :: C4_BETA = 0.1702929d0

    a(:) = C4_BETA
    b(:) = C4_ALPHA
    c(:) = C_1_R
    d(:) = C4_ALPHA
    e(:) = C4_BETA

    return
end subroutine FLT_C4P_CUTOFF_LHS

subroutine FLT_C4P_CUTOFF_RHS(imax, jkmax, u, rhs)
    implicit none

    TINTEGER, intent(in) :: imax, jkmax
    TREAL, intent(in) :: u(jkmax, imax)
    TREAL, intent(out) :: rhs(jkmax, imax)

    integer i
    integer im3, im2, im1, ip1, ip2, ip3

    TREAL, parameter :: C4_A = 0.9891856d0
    TREAL, parameter :: C4_BD2 = 0.66059d0    ! 1.321180 /2
    TREAL, parameter :: C4_CD2 = 0.1666774d0  ! 0.3333548 /2
    TREAL, parameter :: C4_DD2 = 0.679925d-3  ! 0.001359850 /2

    do i = 1, imax
        im3 = MOD(i + imax - 4, imax) + 1   ! i-3
        im2 = MOD(i + imax - 3, imax) + 1   ! i-2
        im1 = MOD(i + imax - 2, imax) + 1   ! i-1
        ip1 = MOD(i, imax) + 1              ! i+1
        ip2 = MOD(i + 1, imax) + 1          ! i+2
        ip3 = MOD(i + 2, imax) + 1          ! i+3
        rhs(:, i) = C4_BD2*(u(:, ip1) + u(:, im1)) + &
                    C4_CD2*(u(:, ip2) + u(:, im2)) + &
                    C4_DD2*(u(:, ip3) + u(:, im3)) + &
                    C4_A*u(:, i)
    end do

    return
end subroutine FLT_C4P_CUTOFF_RHS
