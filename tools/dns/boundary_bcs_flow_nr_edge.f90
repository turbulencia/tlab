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
SUBROUTINE BOUNDARY_BCS_FLOW_NR_EDGE(iflag, jmax, kmax, beta, &
     r, un, v1, v2, p, gama, m1, m2, m3 ,m4, m5, hr, hun, hv1, hv2, he)

  IMPLICIT NONE

#include "types.h"

  TINTEGER iflag, jmax, kmax
  TREAL beta
  
  TREAL, DIMENSION(jmax,*) :: r, un, v1, v2, p, gama
  TREAL, DIMENSION(jmax,*) :: m1, m2, m3, m4, m5
  TREAL, DIMENSION(jmax,*) :: hr, hun, hv1, hv2, he
  
! -------------------------------------------------------------------
  TINTEGER j, k
  TREAL c, Mn, dummy, F1, F2, F3, F4, F5
  
! ###################################################################
! BCs at x_min
! ###################################################################
  IF ( iflag .LT. 0 ) THEN

     DO k = 1,kmax
        DO j = 1,jmax, jmax-1
           c   = SQRT(gama(j,k)*p(j,k)/r(j,k))
           Mn  = un(j,k)/c

           IF ( un(j,k) + c .GT. C_0_R ) THEN
           
! -------------------------------------------------------------------
! Inflow in Ox 
! -------------------------------------------------------------------
              IF ( un(j,k) .GT. C_0_R ) THEN 
                 IF ( j .EQ. 1 ) THEN
                    IF ( v1(j,k) .LT. C_0_R ) THEN ! & Outflow in Oy
                       F1 = C_05_R*m5(j,k)
                       F2 = C_0_R
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = C_0_R
                    ELSE                           ! & Inflow in Oy
                       F1 = C_05_R*m5(j,k)-r(j,k)*c*m2(j,k)
                       F2 = C_0_R
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = C_0_R
                    ENDIF
                 ENDIF
                 IF ( j .EQ. jmax ) THEN 
                    IF ( v1(j,k) .GT. C_0_R ) THEN ! & Outflow in Oy
                       F1 = C_05_R*m1(j,k)
                       F2 = C_0_R
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = C_0_R
                    ELSE                           ! & Inflow in Oy
                       F1 = C_05_R*m1(j,k)-r(j,k)*c*m2(j,k)
                       F2 = C_0_R
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = C_0_R
                    ENDIF
                 ENDIF
                 dummy = ( F2 + C_05_R*(F1+F5) )/c/c

                 hr(j,k)  = hr(j,k)  + dummy
                 hun(j,k) = hun(j,k) + dummy*un(j,k) + ( F5-F1 )*C_05_R/c 
                 hv1(j,k) = hv1(j,k) + dummy*v1(j,k) + r(j,k)*F3
                 hv2(j,k) = hv2(j,k) + dummy*v2(j,k) + r(j,k)*F4
                 he(j,k)  = he(j,k)  + C_05_R*( F1+F5 )/(gama(j,k)-C_1_R)

! -------------------------------------------------------------------
! Outflow in Ox 
! -------------------------------------------------------------------
              ELSE
                 IF ( j .EQ. 1 ) THEN
                    IF ( v1(j,k) .LT. C_0_R ) THEN ! & Outflow in Oy
                       F1 = C_05_R*m5(j,k)
                       F2 = C_0_R
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = beta*C_05_R*m5(j,k)
                    ELSE                           ! & Inflow in Oy
                       F1 = C_05_R*m5(j,k) - r(j,k)*c*m2(j,k)
                       F2 = m3(j,k)
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F4 = m4(j,k)
                       F5 = beta*( C_05_R*m5(j,k) + r(j,k)*c*m2(j,k) )
                    ENDIF
                 ENDIF
                 
                 IF ( j .EQ. jmax ) THEN
                    IF ( v1(j,k) .GT. C_0_R ) THEN ! & Outflow in Oy
                       F1 = C_05_R*m1(j,k)
                       F2 = C_0_R
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = beta*C_05_R*m1(j,k)
                    ELSE                           ! & Inflow in Oy
                       F1 = C_05_R*m1(j,k) - r(j,k)*c*m2(j,k)
                       F2 = m3(j,k)
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F4 = m4(j,k)
                       F5 = beta*( C_05_R*m1(j,k) + r(j,k)*c*m2(j,k) )
                    ENDIF
                 ENDIF
                 dummy = ( F2 + C_05_R*(F1+F5) )/c/c

                 hr(j,k)  = hr(j,k)  + dummy
                 hun(j,k) = hun(j,k) + dummy*un(j,k) + ( F5-F1 )*C_05_R/c 
                 hv1(j,k) = hv1(j,k) + dummy*v1(j,k) + r(j,k)*F3
                 hv2(j,k) = hv2(j,k) + dummy*v2(j,k) + r(j,k)*F4
                 he(j,k)  = he(j,k)  + C_05_R*( F1+F5 )/(gama(j,k)-C_1_R)
                 
              ENDIF
           
           ENDIF ! subsonic branch

        ENDDO
     ENDDO
     
! ###################################################################
! BCs at x_max
! ###################################################################
  ELSE IF ( iflag .GT. 0 ) THEN
     
     DO k = 1,kmax
        DO j = 1,jmax, jmax-1
           c   = SQRT(gama(j,k)*p(j,k)/r(j,k))
           Mn  = un(j,k)/c
           
           IF ( un(j,k) - c .LT. C_0_R ) THEN
           
! -------------------------------------------------------------------
! Inflow in Ox 
! -------------------------------------------------------------------
              IF ( un(j,k) .LT. C_0_R ) THEN
                 dummy = C_0_R

                 IF ( j .EQ. 1 ) THEN
                    IF ( v1(j,k) .LT. C_0_R ) THEN ! & Outflow in Oy
                       F1 = C_0_R
                       F2 = C_0_R
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F3 = C_0_R
                       F4 = C_0_R
                       F5 = C_05_R*m5(j,k)
                    ELSE                           ! & Inflow in Oy
                       F1 = C_0_R
                       F2 = C_0_R
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F3 = C_0_R
                       F4 = C_0_R
                       F5 = C_05_R*m5(j,k)+r(j,k)*c*m2(j,k)
                    ENDIF
                 ENDIF
                 IF ( j .EQ. jmax ) THEN
                    IF ( v1(j,k) .GT. C_0_R ) THEN ! & Outflow in Oy
                       F1 = C_0_R
                       F2 = C_0_R
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F3 = C_0_R
                       F4 = C_0_R
                       F5 = C_05_R*m1(j,k)
                    ELSE                           ! & Inflow in Oy
                       F1 = C_0_R
                       F2 = C_0_R
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F3 = C_0_R
                       F4 = C_0_R
                       F5 = C_05_R*m1(j,k)+r(j,k)*c*m2(j,k)
                    ENDIF
                 ENDIF
                 dummy = ( F2 + C_05_R*(F1+F5) )/c/c

                 hr(j,k)  = hr(j,k)  + dummy
                 hun(j,k) = hun(j,k) + dummy*un(j,k) + ( F5-F1 )*C_05_R/c 
                 hv1(j,k) = hv1(j,k) + dummy*v1(j,k) + r(j,k)*F3
                 hv2(j,k) = hv2(j,k) + dummy*v2(j,k) + r(j,k)*F4
                 he(j,k)  = he(j,k)  + C_05_R*( F1+F5 )/(gama(j,k)-C_1_R)

! -------------------------------------------------------------------
! Outflow in Ox 
! -------------------------------------------------------------------
              ELSE
                 IF ( j .EQ. 1 ) THEN
                    IF ( v1(j,k) .GT. C_0_R ) THEN ! Inflow in Oy
                       F1 = beta*( C_05_R*m5(j,k) - r(j,k)*c*m2(j,k) )
                       F2 = m3(j,k)
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F4 = m4(j,k)
                       F5 = C_05_R*m5(j,k) + r(j,k)*c*m2(j,k)
                    ELSE                           ! Outflow in Oy
                       F1 = beta*C_05_R*m5(j,k)
                       F2 = C_0_R
                       F3 = C_05_R*m5(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = C_05_R*m5(j,k)
                    ENDIF
                 ENDIF
                    
                 IF ( j .EQ. jmax ) THEN
                    IF ( v1(j,k) .LT. C_0_R ) THEN ! Inflow in Oy
                       F1 = beta*( C_05_R*m1(j,k) - r(j,k)*c*m2(j,k) )
                       F2 = m3(j,k)
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F4 = m4(j,k)
                       F5 = C_05_R*m1(j,k) + r(j,k)*c*m2(j,k)
                    ELSE                           ! Outflow in Oy
                       F1 = beta*C_05_R*m1(j,k)
                       F2 = C_0_R
                       F3 =-C_05_R*m1(j,k)/r(j,k)/c
                       F4 = C_0_R
                       F5 = C_05_R*m1(j,k)
                    ENDIF
                 ENDIF
                 dummy = ( F2 + C_05_R*(F1+F5) )/c/c

                 hr(j,k)  = hr(j,k)  + dummy
                 hun(j,k) = hun(j,k) + dummy*un(j,k) + ( F5-F1 )*C_05_R/c 
                 hv1(j,k) = hv1(j,k) + dummy*v1(j,k) + r(j,k)*F3
                 hv2(j,k) = hv2(j,k) + dummy*v2(j,k) + r(j,k)*F4
                 he(j,k)  = he(j,k)  + C_05_R*( F1+F5 )/(gama(j,k)-C_1_R)

              ENDIF
           
           ENDIF ! subsonic branch
        
        ENDDO
     ENDDO
     
  ENDIF
  
  RETURN
END SUBROUTINE BOUNDARY_BCS_FLOW_NR_EDGE
