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
!# transport equations of species mass fractions.
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag   In   Flag selecting BC at xmin (0) or at xmax (1)
!#
!# un      In   Velocity normal to the boundary
!# gn      In   Constant body force normal to the boundary
!#
!########################################################################
      SUBROUTINE BOUNDARY_BCS_SCAL_NR&
           (iflag, nt, pl_const, pl_pref,&
           r, un, z1, p, gama, drdn, dundn, dz1dn, dpdn, gn, hz1)

      IMPLICIT NONE

#include "types.h"

      TINTEGER iflag, nt
      TREAL pl_const, pl_pref
      TREAL r(*), un(*), z1(*), p(*), gama(*)
      TREAL drdn(*), dundn(*), dz1dn(*), dpdn(*), gn
      TREAL hz1(*)

! -------------------------------------------------------------------
      TINTEGER i
      TREAL c, Mn, dummy

! ###################################################################
! BCs at x_min
! ###################################################################
      IF ( iflag .EQ. 0 ) THEN

         DO i = 1, nt            
            c  = SQRT(gama(i)*p(i)/r(i))
            Mn = un(i)/c
            
            IF ( un(i) + c .GT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
               IF ( un(i) .GT. C_0_R ) THEN
                  dummy = C_05_R*( r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R-Mn)/c*dpdn(i) &
                       - r(i)*gn/c )

                  hz1(i) = un(i)*z1(i)*drdn(i) + r(i)*un(i)*dz1dn(i) + dummy*z1(i)
                   
! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
               ELSE
                  dummy = C_05_R*( r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R+Mn)/c*dpdn(i) &
                       - r(i)*gn/c - pl_const*(p(i)-pl_pref)/c )

                  hz1(i) = dummy*z1(i)
                  
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

            IF ( un(i) - c .LT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
               IF ( un(i) .LT. C_0_R ) THEN
                  dummy = C_05_R*( r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R+Mn)/c*dpdn(i) &
                       + r(i)*gn/c )

                  hz1(i) = un(i)*z1(i)*drdn(i) + r(i)*un(i)*dz1dn(i) + dummy*z1(i)
                  
! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
               ELSE
                  dummy = C_05_R*( r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R-Mn)/c*dpdn(i) &
                       + r(i)*gn/c - pl_const*(p(i)-pl_pref)/c)

                  hz1(i) = dummy*z1(i)

               ENDIF
               
            ENDIF ! subsonic branch

         ENDDO
         
      ENDIF

      RETURN
      END
