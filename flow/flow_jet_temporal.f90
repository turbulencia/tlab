!########################################################################
!# Tool/Library FLOW
!#
!########################################################################
!# HISTORY
!#
!# 2007/10/28 - J.P. Mellado
!#              Created. Copied from FLOW_SHEAR_TEMPORAL
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag   In    type of profile:
!#               0 constant value
!#               1 linear
!#               2 hyperbolic tangent
!#               3 error function
!#               4 two linear-varying layers with error function separation
!#
!########################################################################
#include "types.h"

FUNCTION FLOW_JET_TEMPORAL(iflag, thick, delta, mean, diam, ycenter, param, y)

  IMPLICIT NONE

  TINTEGER iflag
  TREAL thick, delta, mean, diam, ycenter, param(*)
  TREAL y, FLOW_JET_TEMPORAL, ERF

! -------------------------------------------------------------------
  TREAL xi, amplify

! ###################################################################
! position relative to yref
  xi = y - ycenter

! symmetry
  IF ( xi .GT. C_0_R ) THEN; xi =  y - ycenter - C_05_R*diam
  ELSE;                      xi =-(y - ycenter + C_05_R*diam); ENDIF

! -------------------------------------------------------------------
! base state varying between two constant levels
! -------------------------------------------------------------------
  IF ( thick .EQ. C_0_R ) THEN
     amplify = C_0_R
     IF      ( xi .LE. C_0_R ) THEN; amplify = C_05_R
     ELSE IF ( xi .GT. C_0_R ) THEN; amplify =-C_05_R; ENDIF

  ELSE
     IF      ( iflag .EQ. 0 ) THEN; amplify = C_0_R
     ELSE IF ( iflag .EQ. 1 ) THEN; amplify =-xi/thick
     ELSE IF ( iflag .EQ. 2 ) THEN; amplify = C_05_R*TANH(-C_05_R*xi/thick)
     ELSE IF ( iflag .EQ. 3 ) THEN; amplify = C_05_R* ERF(-C_05_R*xi/thick)
     ELSE IF ( iflag .EQ. 4 ) THEN; amplify = C_05_R* ERF(-C_05_R*xi/thick); ENDIF
  ENDIF

! mean profile
  FLOW_JET_TEMPORAL = mean + delta*amplify

! -------------------------------------------------------------------
! special profiles
! -------------------------------------------------------------------
! two linear-varying layers
  IF ( iflag .EQ. 4 ) THEN
     IF ( xi .LT. C_0_R ) THEN; FLOW_JET_TEMPORAL = FLOW_JET_TEMPORAL + param(1)*xi
     ELSE;                      FLOW_JET_TEMPORAL = FLOW_JET_TEMPORAL + param(2)*xi; ENDIF
  ENDIF

  RETURN
END FUNCTION FLOW_JET_TEMPORAL
