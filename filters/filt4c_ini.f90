!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2009/01/14 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# 4th-Order Compact Filter from Lele [J. Comp. Phys., V103, 1992]
!#
!# uf_i + alpha*(uf_i-1 + uf_i+1) = a*u_i + b*(u_i-1 + u_i+1) + c*(u_i-2 + u_i+2)
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE FILT4C_INI(imax, i1bc, alpha, dx, cxi)
  
  IMPLICIT NONE
  
#include "types.h"
  
  TINTEGER imax, i1bc

  TREAL cxi(imax,*)
  TREAL dx(imax)
  TREAL alpha

! -----------------------------------------------------------------------
  TREAL ac_loc
  TINTEGER i

! #######################################################################
  IF ( imax .GT. 1 ) THEN
     
! #######################################################################
! Calculate constants for LHS
! #######################################################################
     DO i = 1,imax
        cxi(i,6) = alpha
     ENDDO

! #######################################################################
! Calculate constants for RHS, interior point
! #######################################################################
     ac_loc = (C_5_R + C_6_R*alpha)/C_8_R

     i=1
     cxi(i,1) = (ac_loc-C_1_R)*dx(i)*dx(i+1)*(dx(i+1)+dx(i+2))/&
          (dx(i)*(dx(i) + dx(i)+dx(i+1))*&
          (dx(i)+dx(i)+dx(i+1)+dx(i+2)))
     cxi(i,2) = alpha &
          + (C_1_R-ac_loc)*(dx(i)+dx(1))**2*dx(i+1)*(dx(i+1)+dx(i+2))/&
          (dx(i)*dx(i)*(dx(i)+dx(i+1))*(dx(i)+dx(i)+dx(i+1)+dx(i+2)))&
          + (ac_loc-C_1_R)*(dx(i) + dx(i))*dx(i+1)*(dx(i+1) + dx(i+2))**2/&
          (dx(i)*(dx(i)+dx(i+1))*(dx(i)+dx(i+1)+dx(i+2))*&
          (dx(i) + dx(i) + dx(i+1) + dx(i+2)))
     cxi(i,3) = ac_loc 
     cxi(i,4) = alpha + (C_1_R-ac_loc)*dx(i)&
          *(dx(i)+dx(1))*(dx(i+1)+dx(i+2))/&
          (dx(i+2)*(dx(i)+dx(i+1))*(dx(i)+dx(1)+dx(i+1)))
     cxi(i,5) = (ac_loc-C_1_R)*dx(i)*(dx(i)+dx(i))*dx(i+1)/&
          (dx(i+2)*(dx(i)+dx(i+1)+dx(i+2))*(dx(i)+dx(i)+&
          dx(i+1)+dx(i+2)))

     DO i = 2,imax-2
        cxi(i,1) = (ac_loc-C_1_R)*dx(i)*dx(i+1)*(dx(i+1)+dx(i+2))/&
             (dx(i-1)*(dx(i) + dx(i-1)+dx(i+1))*&
             (dx(i)+dx(i-1)+dx(i+1)+dx(i+2)))
        cxi(i,2) = alpha &
             + (C_1_R-ac_loc)*(dx(i)+dx(i-1))**2*dx(i+1)*(dx(i+1)+dx(i+2))/&
             (dx(i)*dx(i-1)*(dx(i)+dx(i+1))*(dx(i)+dx(i-1)+dx(i+1)+dx(i+2)))&
             + (ac_loc-C_1_R)*(dx(i) + dx(i-1))*dx(i+1)*(dx(i+1)+dx(i+2))**2/&
             (dx(i)*(dx(i)+dx(i+1))*(dx(i)+dx(i+1)+dx(i+2))*&
             (dx(i) + dx(i-1) + dx(i+1) + dx(i+2)))
        cxi(i,3) = ac_loc 
        cxi(i,4) = alpha + (C_1_R-ac_loc)&
             *dx(i)*(dx(i)+dx(i-1))*(dx(i+1)+dx(i+2))/&
             (dx(i+2)*(dx(i)+dx(i+1))*(dx(i)+dx(i-1)+dx(i+1)))
        cxi(i,5) = (ac_loc-C_1_R)*dx(i)*(dx(i)+dx(i-1))*dx(i+1)/&
             (dx(i+2)*(dx(i)+dx(i+1)+dx(i+2))*(dx(i)+dx(i-1)+&
             dx(i+1)+dx(i+2)))
     ENDDO

     i=imax-1
     cxi(i,1) = (ac_loc-C_1_R)*dx(i)*dx(i+1)*(dx(i+1)+dx(i+1))/&
          (dx(i-1)*(dx(i) + dx(i-1)+dx(i+1))*&
          (dx(i)+dx(i-1)+dx(i+1)+dx(i+1)))
     cxi(i,2) = alpha &
          + (C_1_R-ac_loc)*(dx(i)+dx(i-1))**2*dx(i+1)*(dx(i+1)+dx(i+1))/&
          (dx(i)*dx(i-1)*(dx(i)+dx(i+1))*(dx(i)+dx(i-1)+dx(i+1)+dx(i+1)))&
          + (ac_loc-C_1_R)*(dx(i) + dx(i-1))*dx(i+1)*(dx(i+1)+dx(i+1))**2/&
          (dx(i)*(dx(i)+dx(i+1))*(dx(i)+dx(i+1)+dx(i+1))*&
          (dx(i) + dx(i-1) + dx(i+1) + dx(i+1)))
     cxi(i,3) = ac_loc 
     cxi(i,4) = alpha + (C_1_R-ac_loc)&
          *dx(i)*(dx(i)+dx(i-1))*(dx(i+1)+dx(i+1))/&
          (dx(i+1)*(dx(i)+dx(i+1))*(dx(i)+dx(i-1)+dx(i+1)))
     cxi(i,5) = (ac_loc-C_1_R)*dx(i)*(dx(i)+dx(i-1))*dx(i+1)/&
          (dx(i+1)*(dx(i)+dx(i+1)+dx(i+1))*(dx(i)+dx(i-1)+&
          dx(i+1)+dx(i+1)))

     i=imax
     cxi(i,1) = (ac_loc-C_1_R)*dx(i)*dx(i)*(dx(i)+dx(i))/&
          (dx(i-1)*(dx(i) + dx(i-1)+dx(i))*&
          (dx(i)+dx(i-1)+dx(i)+dx(i)))
     cxi(i,2) = alpha &
          + (C_1_R-ac_loc)*(dx(i)+dx(i-1))**2*dx(i)*(dx(i)+dx(i))/&
          (dx(i)*dx(i-1)*(dx(i)+dx(i))*(dx(i)+dx(i-1)+dx(i)+dx(i)))&
          + (ac_loc-C_1_R)*(dx(i) + dx(i-1))*dx(i)*(dx(i) + dx(i))**2/&
          (dx(i)*(dx(i)+dx(i))*(dx(i)+dx(i)+dx(i))*&
          (dx(i) + dx(i-1) + dx(i) + dx(i)))
     cxi(i,3) = ac_loc 
     cxi(i,4) = alpha + (C_1_R-ac_loc)*dx(i)*(dx(i)+dx(i-1))*(dx(i)+dx(i))/&
          (dx(i)*(dx(i)+dx(i))*(dx(i)+dx(i-1)+dx(i)))
     cxi(i,5) = (ac_loc-C_1_R)*dx(i)*(dx(i)+dx(i-1))*dx(i)/&
          (dx(i)*(dx(i)+dx(i)+dx(i))*(dx(i)+dx(i-1)+&
          dx(i)+dx(i)))

! #######################################################################
! biased formulation at the BCs
! #######################################################################
     IF ( i1bc .GT. 0 ) THEN

! -----------------------------------------------------------------------
! i = 1
! -----------------------------------------------------------------------
        i = 1
        ac_loc = (C_15_R + alpha)/C_16_R

        cxi(i,1) = ac_loc
        cxi(i,2)= alpha + (C_1_R-ac_loc)&
             *(dx(2)+dx(3))*(dx(2)+dx(3)+dx(4))*&
             (dx(2)+dx(3)+dx(4)+dx(5))/(dx(3)*(dx(3)+dx(4))*(dx(3)+dx(4)+dx(5)))
        cxi(i,3)= (ac_loc-C_1_R)&
             *dx(2)*(dx(2)+dx(3)+dx(4))*(dx(2)+dx(3)+dx(4)+dx(5))/&
             (dx(3)*dx(4)*(dx(4)+dx(5)))
        cxi(i,4)= (C_1_R-ac_loc)*dx(2)*(dx(2)+dx(3))*(dx(2)+dx(3)+dx(4)+dx(5))/&
             (dx(4)*dx(5)*(dx(3)+dx(4)))
        cxi(i,5)= (ac_loc-C_1_R)*dx(2)*(dx(2)+dx(3))*(dx(2)+dx(3)+dx(4))/&
             (dx(5)*(dx(4)+dx(5))*(dx(3)+dx(4)+dx(5)))

! -----------------------------------------------------------------------
! i = 2
! -----------------------------------------------------------------------
        i = 2
        ac_loc = (C_3_R + C_2_R*alpha)/C_4_R

        cxi(i,1) = alpha + (C_1_R &
             - ac_loc)*dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/&
             ((dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))*dx(5))&
             + (-C_1_R + ac_loc)&
             *dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/&
             ((dx(2) + dx(3))*dx(5)*(dx(2) + dx(3) + dx(4) + dx(5)))
        cxi(i,2) = ac_loc
        cxi(i,3) = alpha + dx(2)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/&
             ((dx(2) + dx(3))*dx(4)*(dx(4) + dx(5)))&
             - ac_loc*dx(2)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/&
             ((dx(2) + dx(3))*dx(4)*(dx(4) + dx(5)))
        cxi(i,4) = (-C_1_R + ac_loc)*dx(2)*dx(3)*(dx(3) + dx(4) + dx(5))/&
             (dx(4)*(dx(2) + dx(3) + dx(4))*dx(5))
        cxi(i,5) = (C_1_R - ac_loc)*dx(2)*dx(3)*(dx(3) + dx(4))/&
             (dx(5)*(dx(4) + dx(5))*(dx(2) + dx(3) + dx(4) + dx(5)))

! -----------------------------------------------------------------------
! i = imax-1
! -----------------------------------------------------------------------
        i = imax-1
        ac_loc = (C_3_R + C_2_R*alpha)/C_4_R

        cxi(i,5) = alpha + (C_1_R - ac_loc)*&
             dx(i)*(dx(i) + dx(i-1))*(dx(i) + dx(i-1) + dx(i-2))/&
             ((dx(i) + dx(i+1))*(dx(i) + dx(i-1) + dx(i+1))*&
             (dx(i) + dx(i-1) + dx(i-2) + dx(i+1)))
        cxi(i,4) = ac_loc
        cxi(i,3) = alpha + (C_1_R - ac_loc)*&
             (dx(i) + dx(i-1))*(dx(i) + dx(i-1) + dx(i-2))*dx(i+1)/&
             (dx(i-1)*dx(i-2)*(dx(i) + dx(i+1))) &
             + (-C_1_R + ac_loc)*(dx(i) + dx(i-1))*&
             (dx(i) + dx(i-1) + dx(i-2))*dx(i+1)/&
             (dx(i-2)*(dx(i-1) + dx(i-2))*(dx(i) + dx(i+1)))
        cxi(i,2) = (-C_1_R + ac_loc)&
             *dx(i)*(dx(i) + dx(i-1) + dx(i-2))*dx(i+1)/&
             (dx(i-1)*dx(i-2)*(dx(i) + dx(i-1) + dx(i+1)))
        cxi(i,1) = (C_1_R - ac_loc)*dx(i)*(dx(i) + dx(i-1))*dx(i+1)/&
             (dx(i-2)*(dx(i-1) + dx(i-2))*(dx(i) + dx(i-1) + dx(i-2) + &
             dx(i+1)))

! -----------------------------------------------------------------------
! i = imax
! -----------------------------------------------------------------------
        i = imax
        ac_loc = (C_15_R + alpha)/C_16_R

        cxi(i,5) = ac_loc
        cxi(i,4) = alpha + (C_1_R-ac_loc)&
             *(dx(i)+dx(i-1))*(dx(i)+dx(i-1)+dx(i-2))*&
             (dx(i)+dx(i-1)+dx(i-2)+dx(i-3))/&
             (dx(i-1)*(dx(i-1)+dx(i-2))*(dx(i-1)+dx(i-2)+dx(i-3)))
        cxi(i,3) = (ac_loc-C_1_R)*dx(i)*(dx(i)+dx(i-1)+dx(i-2))*&
             (dx(i)+dx(i-1)+dx(i-2)+dx(i-3))/&
             (dx(i-1)*dx(i-2)*(dx(i-2)+dx(i-3)))
        cxi(i,2) = (C_1_R-ac_loc)*dx(i)&
             *(dx(i)+dx(i-1))*(dx(i)+dx(i-1)+dx(i-2)+dx(i-3))/&
             (dx(i-2)*dx(i-3)*(dx(i-1)+dx(i-2)))
        cxi(i,1) = (ac_loc-C_1_R)*dx(i)&
             *(dx(i)+dx(i-1))*(dx(i)+dx(i-1)+dx(i-2))/&
             (dx(i-3)*(dx(i-2)+dx(i-3))*(dx(i-1)+dx(i-2)+dx(i-3)))

     ENDIF


  ENDIF

  RETURN
END SUBROUTINE FILT4C_INI

