!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/10/23 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the nonreflective boundary conditions for the 
!# transport equations of density, momentum and internal energy.
!#
!# Poinsot&Lele correction term is included, with pl_out.
!#
!# Copied from BOUNDARY_BCS_FLOW_NR to use in the formulations using the
!# internal energy instead of the total energy.
!#
!# The case of internal energy is like that of p, devided by (\gamma-1), 
!# and with the formation part coming from the scalar equations.
!#
!# Copied from BOUNDARY_BCS_FLOW_NR_1 adding forcing terms for inflow
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag     In   Flag selecting BC at xmin (<0) or at xmax (>0)
!#                1. nonreflective
!#                2. fluctuation
!#                3. mean
!#                4. fluctuation+mean
!# idir      In   Flag to predefined models for relaxation terms
!#                1. OX direction
!#                2. OY direction
!#                3. all terms
!# un        In   Velocity normal to the boundary
!# v1        In   Velocity transversal to the normal
!# v2        In   Velocity transversal to the normal
!# gn        In   Constant body force normal to the boundary
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_FLOW_NR_3(iflag, idir, nt, pl_out, pl_inf, inf_rhs, bf, bf_shape, &
     r, un, v1, v2, p, gama, drdn, dundn, dv1dn, dv2dn, dpdn, gn, hr, hun, hv1, hv2, he)
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER iflag, idir, nt
  TREAL pl_out, pl_inf
  TREAL inf_rhs(nt,*), bf(nt,*), bf_shape(*)
  TREAL r(*), un(*), v1(*), v2(*), p(*), gama(*)
  TREAL drdn(*), dundn(*), dv1dn(*), dv2dn(*), dpdn(*), gn
  TREAL hr(*), hun(*), hv1(*), hv2(*), he(*)

! -------------------------------------------------------------------
  TINTEGER i
  TREAL c, Mn, M2, dummy
  TREAL F1, F2, F3, F4, F5

! ###################################################################
! BCs at x_min
! ###################################################################
  IF ( iflag .LT. 0 ) THEN

     DO i = 1, nt            
        c   = SQRT(gama(i)*p(i)/r(i))
        Mn  = un(i)/c
        M2  = C_05_R*( un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i) )/c/c

        IF ( un(i) + c .GT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .GT. C_0_R ) THEN
              dummy = C_05_R*( r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R-Mn)/c*dpdn(i) - r(i)*gn/c )

              hr(i)  = un(i)*drdn(i) + dummy
              hun(i) = un(i)*un(i)*drdn(i) + dummy*c*(1+Mn) + Mn*dpdn(i)
              hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
              hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
              he(i)  = (un(i)*dpdn(i) + dummy*c*c)/(gama(i)-C_1_R)

! add forcing terms
              F2 = C_0_R
              F3 = C_0_R
              F4 = C_0_R
              F5 = C_0_R
              IF ( ABS(iflag) .EQ. 2 .OR. ABS(iflag) .EQ. 4 ) THEN ! fluctuation
                 F2 = F2 + (inf_rhs(i,1) - inf_rhs(i,5)/c/c    )*bf_shape(i)
                 F3 = F3 + (inf_rhs(i,3)                       )*bf_shape(i)
                 F4 = F4 + (inf_rhs(i,4)                       )*bf_shape(i)
                 F5 = F5 + (inf_rhs(i,5) + inf_rhs(i,2)*r(i)*c )*bf_shape(i)
              ENDIF
              IF ( ABS(iflag) .EQ. 3 .OR. ABS(iflag) .EQ. 4 ) THEN ! mean
                 IF ( idir .EQ. 1 ) THEN      ! OX direction; possible v1 forcing w/ F3
!                    F2 = F2 - pl_inf*( (r(i)-p(i)/c/c) - (bf(i,1)-bf(i,5)/c/c) )
                    F2 = F2 -pl_inf*  ( r(i)  - bf(i,1) )*bf_shape(i) &
                            -pl_out*c*( r(i)  - bf(i,1) )*(C_1_R-bf_shape(i))
                    F3 = F3 -pl_inf*  ( v1(i) - bf(i,3) )*bf_shape(i)
                    F4 = F4 -pl_inf*  ( v2(i) - bf(i,4) )*bf_shape(i) &
                            -pl_out*c*( v2(i) - bf(i,4) )*(C_1_R-bf_shape(i))
                    F5 = F5 -pl_inf*  ( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*bf_shape(i) &
                            -pl_out*c*( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*(C_1_R-bf_shape(i))
                 ELSE IF ( idir .EQ. 2 ) THEN ! OY direction; no un forcing w/ F5
                    F2 = F2 - pl_inf*c*( r(i)  - bf(i,1) )
                    F3 = F3 - pl_inf*c*( v1(i) - bf(i,3) )
                    F4 = F4 - pl_inf*c*( v2(i) - bf(i,4) )
                    F5 = F5 - pl_inf*c*( p(i)  - bf(i,5) )
                 ENDIF
              ENDIF

              dummy = F2 + C_05_R*F5/c/c

              hr(i)  = hr(i)  + dummy
              hun(i) = hun(i) + un(i)*F2 + C_05_R*(Mn+C_1_R)*F5/c
              hv1(i) = hv1(i) + r(i)*F3 + v1(i)*dummy
              hv2(i) = hv2(i) + r(i)*F4 + v2(i)*dummy
              he(i)  = he(i)  + C_05_R*F5/(gama(i)-C_1_R)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              F5 = C_0_R
! Poinsot & Lele term. Multiplication by c done in the definition of dummy below 
              IF ( idir .EQ. 1 ) THEN      ! OX treatment at xmin; complete forcing
                 F5 =-pl_inf/c*( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*bf_shape(i) &
                     -pl_out*  ( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*(C_1_R-bf_shape(i))
              ELSE IF ( idir .EQ. 2 ) THEN ! OY treatment; no un forcing w/ F5
                 F5 =-pl_out*( p(i) - bf(i,5) )
              ENDIF

              dummy = C_05_R*(r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R+Mn)/c*dpdn(i) - r(i)*gn/c + F5/c)

              hr(i)  = dummy
              hun(i) = dummy*c*(C_1_R+Mn)
              hv1(i) = dummy*v1(i)
              hv2(i) = dummy*v2(i)
              he(i)  = dummy*c*c/(gama(i)-C_1_R)

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
        M2  = C_05_R*( un(i)*un(i) + v1(i)*v1(i) + v2(i)*v2(i) )/c/c

        IF ( un(i) - c .LT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .LT. C_0_R ) THEN
              dummy = C_05_R*( r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R+Mn)/c*dpdn(i) + r(i)*gn/c )

              hr(i)  = un(i)*drdn(i) + dummy
              hun(i) = un(i)*un(i)*drdn(i) - (1-Mn)*c*dummy - Mn*dpdn(i)
              hv1(i) = un(i)*v1(i)*drdn(i) + r(i)*un(i)*dv1dn(i) + dummy*v1(i)
              hv2(i) = un(i)*v2(i)*drdn(i) + r(i)*un(i)*dv2dn(i) + dummy*v2(i)
              he(i)  = (un(i)*dpdn(i) + dummy*c*c)/(gama(i)-C_1_R)

! add forcing terms: only mean
              F1 = C_0_R
              F2 = C_0_R
              F3 = C_0_R
              F4 = C_0_R
              IF ( ABS(iflag) .EQ. 3 .OR. ABS(iflag) .EQ. 4 ) THEN ! mean
                 IF ( idir .EQ. 1 ) THEN      ! OX direction
                    F1 = F1 - pl_inf*c*( (p(i)-r(i)*c*un(i)) - (bf(i,5)-r(i)*c*bf(i,2)) )
                    F2 = F2 - pl_inf*c*( r(i)  - bf(i,1) )
                    F3 = F3 - pl_inf*c*( v1(i) - bf(i,3) )
                    F4 = F4 - pl_inf*c*( v2(i) - bf(i,4) )
                 ELSE IF ( idir .EQ. 2 ) THEN ! OY direction; no un forcing w/ F1
                    F1 = F1 - pl_inf*c*( p(i)  - bf(i,5) )
                    F2 = F2 - pl_inf*c*( r(i)  - bf(i,1) )
                    F3 = F3 - pl_inf*c*( v1(i) - bf(i,3) )
                    F4 = F4 - pl_inf*c*( v2(i) - bf(i,4) )
                 ENDIF
              ENDIF

              dummy = F2 + C_05_R*F1/c/c

              hr(i)  = hr(i)  + dummy
              hun(i) = hun(i) + un(i)*F2 + C_05_R*(Mn-C_1_R)*F1/c
              hv1(i) = hv1(i) + r(i)*F3 + v1(i)*dummy
              hv2(i) = hv2(i) + r(i)*F4 + v2(i)*dummy
              he(i)  = he(i)  + C_05_R*F1/(gama(i)-C_1_R)

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              F1 = C_0_R
! Poinsot & Lele term. Multiplication by c done in the definition of dummy below 
              IF ( idir .EQ. 1 ) THEN      ! OX treatment; no un forcing w/ F1
                 F1 =-pl_out*( p(i)-bf(i,5) )
              ELSE IF ( idir .EQ. 2 ) THEN ! OY treatment; no un forcing w/ F1
                 F1 =-pl_out*( p(i)-bf(i,5) )
              ENDIF
              dummy = C_05_R*(r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R-Mn)/c*dpdn(i) + r(i)*gn/c + F1/c)

              hr(i)  = dummy
              hun(i) =-dummy*c*(C_1_R-Mn)
              hv1(i) = dummy*v1(i)
              hv2(i) = dummy*v2(i)
              he(i)  = dummy*c*c/(gama(i)-C_1_R)

           ENDIF

        ENDIF ! subsonic branch

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BCS_FLOW_NR_3
