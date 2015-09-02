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
SUBROUTINE BOUNDARY_BCS_FLOW_NR_4(iflag, idir, nt, beta, &
     r, un, v1, v2, p, gama, t1, t2, t3, t4, t5, m1, m5, hr, hun, hv1, hv2, he)

  IMPLICIT NONE

#include "types.h"

  TINTEGER iflag, idir, nt
  TREAL beta
  
  TREAL, DIMENSION(*) :: r, un, v1, v2, p, gama
  TREAL, DIMENSION(*) :: t1, t2, t3, t4, t5, m1, m5
  TREAL, DIMENSION(*) :: hr, hun, hv1, hv2, he
  
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

              hr(i)  = hr(i)  + dummy
              hun(i) = hun(i) + C_05_R*(Mn-C_1_R)*t5(i)/c-C_05_R*r(i)*(Mn+C_1_R)*t2(i)-t1(i)*un(i)
              hv1(i) = hv1(i) + dummy*v1(i) - r(i)*t3(i)
              hv2(i) = hv2(i) + dummy*v2(i) - r(i)*t4(i)
              he(i)  = he(i)  - C_05_R*(t5(i)+r(i)*c*t2(i))/(gama(i)-C_1_R)

! recover lateral term for v1 velocity at inflow in Ox
              IF ( idir .EQ. 1 .OR. idir .EQ. 2 ) THEN
                 hv1(i) = hv1(i) - C_05_R*(m5(i)-m1(i))/c
              ENDIF

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              dummy =-C_05_R*(C_1_R-beta)*( r(i)*c*t2(i) + t5(i) )/c/c
              
              hr(i)  = hr(i)  + dummy
              hun(i) = hun(i) + dummy*c*(C_1_R+Mn)
              hv1(i) = hv1(i) + dummy*v1(i)
              hv2(i) = hv2(i) + dummy*v2(i)
              he(i)  = he(i)  + dummy*c*c/(gama(i)-C_1_R)

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

              hr(i)  = hr(i)  + dummy
              hun(i) = hun(i) + C_05_R*(Mn+C_1_R)*t5(i)/c+C_05_R*r(i)*(Mn-C_1_R)*t2(i)-t1(i)*un(i)
              hv1(i) = hv1(i) + dummy*v1(i) - r(i)*t3(i)
              hv2(i) = hv2(i) + dummy*v2(i) - r(i)*t4(i)
              he(i)  = he(i)  - C_05_R*(t5(i)-r(i)*c*t2(i))/(gama(i)-C_1_R)

! recover lateral term for v1 velocity at inflow in Ox
              IF ( idir .EQ. 1 .OR. idir .EQ. 2 ) THEN
                 hv1(i) = hv1(i) - C_05_R*(m5(i)-m1(i))/c
              ENDIF

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              dummy = C_05_R*(C_1_R-beta)*( r(i)*c*t2(i) - t5(i) )/c/c
              
              hr(i)  = hr(i)  + dummy
              hun(i) = hun(i) - dummy*c*(C_1_R-Mn)
              hv1(i) = hv1(i) + dummy*v1(i)
              hv2(i) = hv2(i) + dummy*v2(i)
              he(i)  = he(i)  + dummy*c*c/(gama(i)-C_1_R)
           ENDIF
           
        ENDIF ! subsonic branch
        
     ENDDO
     
  ENDIF
  
  RETURN
END SUBROUTINE BOUNDARY_BCS_FLOW_NR_4
