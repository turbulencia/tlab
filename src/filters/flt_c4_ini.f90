#include "types.h"

!########################################################################
!# DESCRIPTION
!#
!# 4th-Order Compact Filter from Lele [J. Comp. Phys., V103, 1992]
!#
!# uf_i + alpha*(uf_i-1 + uf_i+1) = a*u_i + b*(u_i-1 + u_i+1) + c*(u_i-2 + u_i+2)
!#
!########################################################################
subroutine FLT_C4_INI(dx, f)

    use TLAB_TYPES, only: filter_dt

    implicit none

    TREAL, dimension(*), intent(IN) :: dx
    type(filter_dt), intent(INOUT) :: f

! -----------------------------------------------------------------------
    TREAL ac_loc, alpha
    TINTEGER i

! #######################################################################
    if (f%size > 1) then

        alpha = f%parameters(1)

! #######################################################################
! Calculate constants for LHS
! #######################################################################
        do i = 1, f%size
            f%coeffs(i, 6) = alpha
        end do

! #######################################################################
! Calculate constants for RHS, interior point
! #######################################################################
        ac_loc = (C_5_R + C_6_R*alpha)/C_8_R

        i = 1
        f%coeffs(i, 1) = (ac_loc - C_1_R)*dx(i)*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                         (dx(i)*(dx(i) + dx(i) + dx(i + 1))* &
                          (dx(i) + dx(i) + dx(i + 1) + dx(i + 2)))
        f%coeffs(i, 2) = alpha &
                         + (C_1_R - ac_loc)*(dx(i) + dx(1))**2*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                         (dx(i)*dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i) + dx(i + 1) + dx(i + 2))) &
                         + (ac_loc - C_1_R)*(dx(i) + dx(i))*dx(i + 1)*(dx(i + 1) + dx(i + 2))**2/ &
                         (dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i + 1) + dx(i + 2))* &
                          (dx(i) + dx(i) + dx(i + 1) + dx(i + 2)))
        f%coeffs(i, 3) = ac_loc
        f%coeffs(i, 4) = alpha + (C_1_R - ac_loc)*dx(i) &
                         *(dx(i) + dx(1))*(dx(i + 1) + dx(i + 2))/ &
                         (dx(i + 2)*(dx(i) + dx(i + 1))*(dx(i) + dx(1) + dx(i + 1)))
        f%coeffs(i, 5) = (ac_loc - C_1_R)*dx(i)*(dx(i) + dx(i))*dx(i + 1)/ &
                         (dx(i + 2)*(dx(i) + dx(i + 1) + dx(i + 2))*(dx(i) + dx(i) + &
                                                                     dx(i + 1) + dx(i + 2)))

        do i = 2, f%size - 2
            f%coeffs(i, 1) = (ac_loc - C_1_R)*dx(i)*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                             (dx(i - 1)*(dx(i) + dx(i - 1) + dx(i + 1))* &
                              (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2)))
            f%coeffs(i, 2) = alpha + (C_1_R - ac_loc) &
                             *(dx(i) + dx(i - 1))**2*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                             (dx(i)*dx(i - 1)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2))) &
                             + (ac_loc - C_1_R)*(dx(i) + dx(i - 1))*dx(i + 1)*(dx(i + 1) + dx(i + 2))**2/ &
                             (dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i + 1) + dx(i + 2))* &
                              (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2)))
            f%coeffs(i, 3) = ac_loc
            f%coeffs(i, 4) = alpha + (C_1_R - ac_loc) &
                             *dx(i)*(dx(i) + dx(i - 1))*(dx(i + 1) + dx(i + 2))/ &
                             (dx(i + 2)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1)))
            f%coeffs(i, 5) = (ac_loc - C_1_R)*dx(i)*(dx(i) + dx(i - 1))*dx(i + 1)/ &
                             (dx(i + 2)*(dx(i) + dx(i + 1) + dx(i + 2))* &
                              (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2)))
        end do

        i = f%size - 1
        f%coeffs(i, 1) = (ac_loc - C_1_R)*dx(i)*dx(i + 1)*(dx(i + 1) + dx(i + 1))/ &
                         (dx(i - 1)*(dx(i) + dx(i - 1) + dx(i + 1))* &
                          (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 1)))
        f%coeffs(i, 2) = alpha &
                         + (C_1_R - ac_loc)*(dx(i) + dx(i - 1))**2*dx(i + 1)*(dx(i + 1) + dx(i + 1))/ &
                         (dx(i)*dx(i - 1)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 1))) &
                         + (ac_loc - C_1_R)*(dx(i) + dx(i - 1))*dx(i + 1)*(dx(i + 1) + dx(i + 1))**2/ &
                         (dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i + 1) + dx(i + 1))* &
                          (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 1)))
        f%coeffs(i, 3) = ac_loc
        f%coeffs(i, 4) = alpha + (C_1_R - ac_loc) &
                         *dx(i)*(dx(i) + dx(i - 1))*(dx(i + 1) + dx(i + 1))/ &
                         (dx(i + 1)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1)))
        f%coeffs(i, 5) = (ac_loc - C_1_R)*dx(i)*(dx(i) + dx(i - 1))*dx(i + 1)/ &
                         (dx(i + 1)*(dx(i) + dx(i + 1) + dx(i + 1))*(dx(i) + dx(i - 1) + &
                                                                     dx(i + 1) + dx(i + 1)))

        i = f%size
        f%coeffs(i, 1) = (ac_loc - C_1_R)*dx(i)*dx(i)*(dx(i) + dx(i))/ &
                         (dx(i - 1)*(dx(i) + dx(i - 1) + dx(i))* &
                          (dx(i) + dx(i - 1) + dx(i) + dx(i)))
        f%coeffs(i, 2) = alpha &
                         + (C_1_R - ac_loc)*(dx(i) + dx(i - 1))**2*dx(i)*(dx(i) + dx(i))/ &
                         (dx(i)*dx(i - 1)*(dx(i) + dx(i))*(dx(i) + dx(i - 1) + dx(i) + dx(i))) &
                         + (ac_loc - C_1_R)*(dx(i) + dx(i - 1))*dx(i)*(dx(i) + dx(i))**2/ &
                         (dx(i)*(dx(i) + dx(i))*(dx(i) + dx(i) + dx(i))* &
                          (dx(i) + dx(i - 1) + dx(i) + dx(i)))
        f%coeffs(i, 3) = ac_loc
        f%coeffs(i, 4) = alpha + (C_1_R - ac_loc)*dx(i)*(dx(i) + dx(i - 1))*(dx(i) + dx(i))/ &
                         (dx(i)*(dx(i) + dx(i))*(dx(i) + dx(i - 1) + dx(i)))
        f%coeffs(i, 5) = (ac_loc - C_1_R)*dx(i)*(dx(i) + dx(i - 1))*dx(i)/ &
                         (dx(i)*(dx(i) + dx(i) + dx(i))*(dx(i) + dx(i - 1) + &
                                                         dx(i) + dx(i)))

! #######################################################################
! biased formulation at the BCs
! #######################################################################
        if (.not. f%periodic) then

! -----------------------------------------------------------------------
! i = 1
! -----------------------------------------------------------------------
            i = 1
            ac_loc = (C_15_R + alpha)/C_16_R

            f%coeffs(i, 1) = ac_loc
            f%coeffs(i, 2) = alpha + (C_1_R - ac_loc) &
                             *(dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))* &
                             (dx(2) + dx(3) + dx(4) + dx(5))/(dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5)))
            f%coeffs(i, 3) = (ac_loc - C_1_R) &
                             *dx(2)*(dx(2) + dx(3) + dx(4))*(dx(2) + dx(3) + dx(4) + dx(5))/ &
                             (dx(3)*dx(4)*(dx(4) + dx(5)))
            f%coeffs(i, 4) = (C_1_R - ac_loc)*dx(2)*(dx(2) + dx(3))*(dx(2) + dx(3) + dx(4) + dx(5))/ &
                             (dx(4)*dx(5)*(dx(3) + dx(4)))
            f%coeffs(i, 5) = (ac_loc - C_1_R)*dx(2)*(dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))/ &
                             (dx(5)*(dx(4) + dx(5))*(dx(3) + dx(4) + dx(5)))

! -----------------------------------------------------------------------
! i = 2
! -----------------------------------------------------------------------
            i = 2
            ac_loc = (C_3_R + C_2_R*alpha)/C_4_R

            f%coeffs(i, 1) = alpha + (C_1_R &
                                      - ac_loc)*dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                             ((dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))*dx(5)) &
                             + (-C_1_R + ac_loc) &
                             *dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                             ((dx(2) + dx(3))*dx(5)*(dx(2) + dx(3) + dx(4) + dx(5)))
            f%coeffs(i, 2) = ac_loc
            f%coeffs(i, 3) = alpha + dx(2)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                             ((dx(2) + dx(3))*dx(4)*(dx(4) + dx(5))) &
                             - ac_loc*dx(2)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                             ((dx(2) + dx(3))*dx(4)*(dx(4) + dx(5)))
            f%coeffs(i, 4) = (-C_1_R + ac_loc)*dx(2)*dx(3)*(dx(3) + dx(4) + dx(5))/ &
                             (dx(4)*(dx(2) + dx(3) + dx(4))*dx(5))
            f%coeffs(i, 5) = (C_1_R - ac_loc)*dx(2)*dx(3)*(dx(3) + dx(4))/ &
                             (dx(5)*(dx(4) + dx(5))*(dx(2) + dx(3) + dx(4) + dx(5)))

! -----------------------------------------------------------------------
! i = f%size-1
! -----------------------------------------------------------------------
            i = f%size - 1
            ac_loc = (C_3_R + C_2_R*alpha)/C_4_R

            f%coeffs(i, 5) = alpha + (C_1_R - ac_loc)* &
                             dx(i)*(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))/ &
                             ((dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1))* &
                              (dx(i) + dx(i - 1) + dx(i - 2) + dx(i + 1)))
            f%coeffs(i, 4) = ac_loc
            f%coeffs(i, 3) = alpha + (C_1_R - ac_loc)* &
                             (dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))*dx(i + 1)/ &
                             (dx(i - 1)*dx(i - 2)*(dx(i) + dx(i + 1))) &
                             + (-C_1_R + ac_loc)*(dx(i) + dx(i - 1))* &
                             (dx(i) + dx(i - 1) + dx(i - 2))*dx(i + 1)/ &
                             (dx(i - 2)*(dx(i - 1) + dx(i - 2))*(dx(i) + dx(i + 1)))
            f%coeffs(i, 2) = (-C_1_R + ac_loc) &
                             *dx(i)*(dx(i) + dx(i - 1) + dx(i - 2))*dx(i + 1)/ &
                             (dx(i - 1)*dx(i - 2)*(dx(i) + dx(i - 1) + dx(i + 1)))
            f%coeffs(i, 1) = (C_1_R - ac_loc)*dx(i)*(dx(i) + dx(i - 1))*dx(i + 1)/ &
                             (dx(i - 2)*(dx(i - 1) + dx(i - 2))*(dx(i) + dx(i - 1) + dx(i - 2) + &
                                                                 dx(i + 1)))

! -----------------------------------------------------------------------
! i = f%size
! -----------------------------------------------------------------------
            i = f%size
            ac_loc = (C_15_R + alpha)/C_16_R

            f%coeffs(i, 5) = ac_loc
            f%coeffs(i, 4) = alpha + (C_1_R - ac_loc) &
                             *(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))* &
                             (dx(i) + dx(i - 1) + dx(i - 2) + dx(i - 3))/ &
                             (dx(i - 1)*(dx(i - 1) + dx(i - 2))*(dx(i - 1) + dx(i - 2) + dx(i - 3)))
            f%coeffs(i, 3) = (ac_loc - C_1_R)*dx(i)*(dx(i) + dx(i - 1) + dx(i - 2))* &
                             (dx(i) + dx(i - 1) + dx(i - 2) + dx(i - 3))/ &
                             (dx(i - 1)*dx(i - 2)*(dx(i - 2) + dx(i - 3)))
            f%coeffs(i, 2) = (C_1_R - ac_loc)*dx(i) &
                             *(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2) + dx(i - 3))/ &
                             (dx(i - 2)*dx(i - 3)*(dx(i - 1) + dx(i - 2)))
            f%coeffs(i, 1) = (ac_loc - C_1_R)*dx(i) &
                             *(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))/ &
                             (dx(i - 3)*(dx(i - 2) + dx(i - 3))*(dx(i - 1) + dx(i - 2) + dx(i - 3)))

        end if

    end if

    return
end subroutine FLT_C4_INI
