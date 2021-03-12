!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2008/11/23 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the transverse terms of nonreflective boundary 
!# conditions for the transport equations of density, momentum and 
!# internal energy.
!#
!# After Lodato et al, JCP 227 (2008), 5105-5143
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_SCAL_NR_4(iflag, nt, beta, r, un, z1, p, gama, t1, t2, t5, tz1, hz1)

  IMPLICIT NONE

#include "types.h"

  TINTEGER iflag, nt
  TREAL beta
  
  TREAL, DIMENSION(*) :: r, un, z1, p, gama
  TREAL, DIMENSION(*) :: t1, t2, t5, tz1
  TREAL, DIMENSION(*) :: hz1
  
! -------------------------------------------------------------------
  TINTEGER i
  TREAL c, Mn, dummy
  
! ###################################################################
! BCs at x_min
! ###################################################################
  IF ( iflag .LT. 0 ) THEN
     
     DO i = 1, nt            
        c   = SQRT(gama(i)*p(i)/r(i))
        Mn  = un(i)/c

        IF ( un(i) + c .GT. C_0_R ) THEN
           
! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .GT. C_0_R ) THEN
              dummy = C_05_R*t5(i)/c/c - C_05_R*r(i)*t2(i)/c - t1(i)

              hz1(i) = hz1(i) + dummy*z1(i) - r(i)*tz1(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              dummy =-C_05_R*(C_1_R-beta)*( r(i)*c*t2(i) + t5(i) )/c/c
              
              hz1(i) = hz1(i) + dummy*z1(i)

           ENDIF
           
        ENDIF ! subsonic branch
        
     ENDDO
     
! ###################################################################
! BCs at x_max
! ###################################################################
  ELSE IF ( iflag .GT. 0 ) THEN
     
     DO i = 1, nt            
        c   = SQRT(gama(i)*p(i)/r(i))
        Mn  = un(i)/c

        IF ( un(i) - c .LT. C_0_R ) THEN
           
! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .LT. C_0_R ) THEN
              dummy = C_05_R*t5(i)/c/c + C_05_R*r(i)*t2(i)/c - t1(i)

              hz1(i) = hz1(i) + dummy*z1(i) - r(i)*tz1(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              dummy = C_05_R*(C_1_R-beta)*( r(i)*c*t2(i) - t5(i) )/c/c
              
              hz1(i) = hz1(i) + dummy*z1(i)

           ENDIF
           
        ENDIF ! subsonic branch
        
     ENDDO
     
  ENDIF
  
  RETURN
END SUBROUTINE BOUNDARY_BCS_SCAL_NR_4
