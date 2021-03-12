!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2003/10/23 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the nonreflective boundary conditions for the 
!# transport equations of species mass fractions.
!#
!# Copied from BOUNDARY_BCS_FLOW_NR_1 adding forcing terms for inflow
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag   In   Flag selecting BC at xmin (<0) or at xmax (>0)
!#              1. nonreflective
!#              2. fluctuation
!#              3. mean
!#              4. fluctuation+mean
!# idir    In   Flag to predefined models, like OX or OY directions
!# un      In   Velocity normal to the boundary
!# gn      In   Constant body force normal to the boundary
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_SCAL_NR_3(iflag, idir, nt, pl_out, pl_inf, inf_rhs, inf_rhs_z, bf, bf_z,&
     bf_shape, r, un, z1, p, gama, drdn, dundn, dz1dn, dpdn, gn, hz1)
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER iflag, idir, nt
  TREAL pl_out, pl_inf
  TREAL inf_rhs(nt,*), inf_rhs_z(nt), bf(nt,*), bf_z(nt), bf_shape(*)
  TREAL r(*), un(*), z1(*), p(*), gama(*)
  TREAL drdn(*), dundn(*), dz1dn(*), dpdn(*), gn
  TREAL hz1(*)

! -------------------------------------------------------------------
  TINTEGER i
  TREAL c, Mn, dummy
  TREAL F1, F2, F5, FZ

! ###################################################################
! BCs at x_min
! ###################################################################
  IF ( iflag .LT. 0 ) THEN

     DO i = 1, nt            
        c  = SQRT(gama(i)*p(i)/r(i))
        Mn = un(i)/c

        IF ( un(i) + c .GT. C_0_R ) THEN

! -------------------------------------------------------------------
! Inflow
! -------------------------------------------------------------------
           IF ( un(i) .GT. C_0_R ) THEN
              dummy = C_05_R*(r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R-Mn)/c*dpdn(i) - r(i)*gn/c)

              hz1(i) = un(i)*z1(i)*drdn(i) + r(i)*un(i)*dz1dn(i) + dummy*z1(i)

! add forcing terms
              F2 = C_0_R
              F5 = C_0_R
              FZ = C_0_R
              IF ( ABS(iflag) .EQ. 2 .OR. ABS(iflag) .EQ. 4 ) THEN ! fluctuation
                 F2 = F2 + (inf_rhs(i,1) - inf_rhs(i,5)/c/c    )*bf_shape(i)
                 F5 = F5 + (inf_rhs(i,5) + inf_rhs(i,2)*r(i)*c )*bf_shape(i)
                 FZ = FZ + (inf_rhs_z(i)                       )*bf_shape(i)
              ENDIF
              IF ( ABS(iflag) .EQ. 3 .OR. ABS(iflag) .EQ. 4 ) THEN ! mean
                 IF ( idir .EQ. 1 ) THEN      ! OX direction
!                    F2 = F2 - pl_inf*( (r(i)-p(i)/c/c) - (bf(i,1)-bf(i,5)/c/c) )
                    F2 = F2 -pl_inf*  ( r(i)  - bf(i,1) )*bf_shape(i) &
                            -pl_out*c*( r(i)  - bf(i,1) )*(C_1_R-bf_shape(i))
                    F5 = F5 -pl_inf*  ( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*bf_shape(i) &
                            -pl_out*c*( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*(C_1_R-bf_shape(i))
                    FZ = FZ -pl_inf*  ( z1(i) - bf_z(i) )*bf_shape(i) &
                            -pl_out*c*( z1(i) - bf_z(i) )*(C_1_R-bf_shape(i))
                 ELSE IF ( idir .EQ. 2 ) THEN ! OY direction; no un forcing w/ F5
                    F2 = F2 - pl_inf*c*( r(i)  - bf(i,1) )
                    F5 = F5 - pl_inf*c*( p(i)  - bf(i,5) )
                    FZ = FZ - pl_inf*c*( z1(i) - bf_z(i) )
                 ENDIF
              ENDIF

              dummy = F2 + C_05_R*F5/c/c

              hz1(i) = hz1(i) + r(i)*FZ + z1(i)*dummy

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              F5 = C_0_R
! Poinsot & Lele term. Multiplication by c is done in the definition of dummy below 
              IF ( idir .EQ. 1 ) THEN      ! OX treatment at xmin; complete forcing
                 F5 =-pl_inf/c*( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*bf_shape(i) &
                     -pl_out*  ( p(i)+r(i)*c*un(i) - (bf(i,5)+r(i)*c*bf(i,2)) )*(C_1_R-bf_shape(i))
              ELSE IF ( idir .EQ. 2 ) THEN ! OY treatment; no un forcing w/ F5
                 F5 =-pl_out*( p(i) - bf(i,5) )
              ENDIF
              dummy = C_05_R*(r(i)*(C_1_R+Mn)*dundn(i) + (C_1_R+Mn)/c*dpdn(i) - r(i)*gn/c + F5/c)

              hz1(i) = dummy*z1(i)

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
              dummy = C_05_R*(r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R+Mn)/c*dpdn(i) + r(i)*gn/c)

              hz1(i) = un(i)*z1(i)*drdn(i) + r(i)*un(i)*dz1dn(i) + dummy*z1(i)

! add forcing terms: only mean
              F1 = C_0_R
              F2 = C_0_R
              FZ = C_0_R
              IF ( ABS(iflag) .EQ. 3 .OR. ABS(iflag) .EQ. 4 ) THEN ! mean
                 IF ( idir .EQ. 1 ) THEN      ! OX direction; only acoustic, with un
                    F1 = F1 - pl_inf*c*( (p(i)-r(i)*c*un(i))-(bf(i,5)-r(i)*c*bf(i,2)))
                    F2 = F2 - pl_inf*c*( r(i)  - bf(i,1) )
                    FZ = FZ - pl_inf*c*( z1(i) - bf_z(i) )
                 ELSE IF ( idir .EQ. 2 ) THEN ! OY direction
                    F1 = F1 - pl_inf*c*( p(i) - bf(i,5) )
                    F2 = F2 - pl_inf*c*( r(i)  - bf(i,1) )
                    FZ = FZ - pl_inf*c*( z1(i) - bf_z(i) )
                 ENDIF
              ENDIF

              dummy = F2 + C_05_R*F1/c/c

              hz1(i) = hz1(i) + r(i)*FZ + z1(i)*dummy

! -------------------------------------------------------------------
! Outflow
! -------------------------------------------------------------------
           ELSE
              F1 = C_0_R
! Poinsot & Lele term. Multiplication by c is done in the definition 
! of dummy below 
              IF ( idir .EQ. 1 ) THEN      ! OX treatment; no un forcing w/ F1
                 F1 =-pl_out*( p(i) - bf(i,5) )
              ELSE IF ( idir .EQ. 2 ) THEN ! OY treatment; no un forcing w/ F1
                 F1 =-pl_out*( p(i) - bf(i,5) )
              ENDIF
              dummy = C_05_R*(r(i)*(C_1_R-Mn)*dundn(i) - (C_1_R-Mn)/c*dpdn(i) + r(i)*gn/c + F1/c)

              hz1(i) = dummy*z1(i)

           ENDIF

        ENDIF ! subsonic branch

     ENDDO

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_BCS_SCAL_NR_3
