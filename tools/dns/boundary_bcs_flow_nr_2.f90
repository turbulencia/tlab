!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2003/06/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the nonreflective boundary conditions for the 
!# transport equations of density, momentum and total energy.
!#
!# Poinsot&Lele correction term is included, with pl_const.
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag   In   Flag selecting BC at xmin (0) or at xmax (1)
!#
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!# gn        In   Constant body force normal to the boundary
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_FLOW_NR_2&
     (iflag, nt, pl_const, pl_pref, &
     r, un, v1, v2, p, gama, drdn, dundn, dv1dn, dv2dn, dpdn, gn,&
     hr, hun, hv1, hv2, he)

  IMPLICIT NONE

#include "types.h"

  TINTEGER iflag, nt
  TREAL pl_const, pl_pref
  TREAL r(*), un(*), v1(*), v2(*), p(*), gama(*)
  TREAL drdn(*), dundn(*), dv1dn(*), dv2dn(*), dpdn(*), gn
  TREAL hr(*), hun(*), hv1(*), hv2(*), he(*)

! -------------------------------------------------------------------
  TINTEGER i
  TREAL c, Mn, M2, dummy

! ###################################################################
! BCs at x_min
! ###################################################################
  IF ( iflag .EQ. 0 ) THEN

     DO i = 1, nt            
        c   = SQRT(gama(i)*p(i)/r(i))
        Mn  = un(i)/c
        M2  = C_05_R*( un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i) )/c/c

        IF ( un(i) + c .GT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .GT. C_0_R ) THEN
              dummy = C_05_R*( r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R-Mn)/c*dpdn(i) &
                   - r(i)*gn/c )

              hr(i)  = un(i)*drdn(i) + dummy
              hun(i) = un(i)*un(i)*drdn(i) + dummy*c*(1+Mn) + Mn*dpdn(i)
              hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
              hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
              he(i)  = un(i)*M2*c*c*drdn(i)+ r(i)*un(i)*(v1(i)*dv1dn(i)+v2(i)*dv2dn(i))&
                   + dummy*c*c*( C_1_R/(gama(i)-C_1_R) + M2 + Mn ) &
                   + un(i)*( C_1_R/(gama(i)-C_1_R) + Mn )*dpdn(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              dummy = C_05_R*( r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R+Mn)/c*dpdn(i) &
                   - r(i)*gn/c - pl_const*(p(i)-pl_pref)/c )

              hr(i)  = dummy
              hun(i) = dummy*c*(C_1_R+Mn)
              hv1(i) = dummy*v1(i)
              hv2(i) = dummy*v2(i)
              he(i)  = dummy*c*c*( C_1_R/(gama(i)-C_1_R) + M2 + Mn )

           ENDIF

        ENDIF ! subsonic branch

     ENDDO

! ###################################################################
! BCs at x_max
! ###################################################################
  ELSE IF ( iflag .EQ. 1 ) THEN

     DO i = 1, nt            
        c   = SQRT(gama(i)*p(i)/r(i))
        Mn  = un(i)/c
        M2  = C_05_R*( un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i) )/c/c

        IF ( un(i) - c .LT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .LT. C_0_R ) THEN
              dummy = C_05_R*( r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R+Mn)/c*dpdn(i) &
                   + r(i)*gn/c )

              hr(i)  = un(i)*drdn(i) + dummy
              hun(i) = un(i)*un(i)*drdn(i) - (1-Mn)*c*dummy - Mn*dpdn(i)
              hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
              hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
              he(i)  = un(i)*M2*c*c*drdn(i)+ r(i)*un(i)*(v1(i)*dv1dn(i)+v2(i)*dv2dn(i))&
                   + dummy*c*c*( C_1_R/(gama(i)-C_1_R) + M2 - Mn ) &
                   + un(i)*( C_1_R/(gama(i)-C_1_R) - Mn )*dpdn(i)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              dummy = C_05_R*( r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R-Mn)/c*dpdn(i) &
                   + r(i)*gn/c - pl_const*(p(i)-pl_pref)/c)

              hr(i)  = dummy
              hun(i) =-dummy*c*(C_1_R-Mn)
              hv1(i) = dummy*v1(i)
              hv2(i) = dummy*v2(i)
              he(i)  = dummy*c*c*( C_1_R/(gama(i)-C_1_R) + M2 - Mn )

           ENDIF

        ENDIF ! subsonic branch

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BCS_FLOW_NR_2
