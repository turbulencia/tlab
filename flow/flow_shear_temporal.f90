#include "types.h"
#include "dns_const.h"
  
!########################################################################
!# Tool/Library FLOW
!#
!########################################################################
!# HISTORY
!#
!# 2007/05/08 - J.P. Mellado
!#              Created
!# 2007/05/30 - J.P. Mellado
!#              Adding case 4. Adding param as argument
!# 2011/01/19 - Cedrick Ansorge
!#              Adding case 7. - Ekman Layer initial perturbation
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
!#               5 Ekman layer U
!#               6 Ekman layer W
!#               7 Ekman layer Perturbed U
!#               8 Bickley jet
!#
!########################################################################
FUNCTION FLOW_SHEAR_TEMPORAL(iflag, thick, delta, mean, ycenter, param, y)

  IMPLICIT NONE
  
  TINTEGER iflag
  TREAL thick, delta, mean, ycenter, param(*)
  TREAL y, FLOW_SHEAR_TEMPORAL
  
! -------------------------------------------------------------------
  TREAL yrel, xi, amplify, zamp, cnought !, ERF

! ###################################################################
  yrel = y - ycenter ! position relative to ycenter
  amplify = C_0_R    ! default

! -------------------------------------------------------------------
! base state varying between two constant levels
! -------------------------------------------------------------------
  IF ( thick .EQ. C_0_R ) THEN
     IF ( iflag .GT. 0 ) THEN 
     IF      ( yrel .LE. C_0_R ) THEN; amplify = C_05_R
     ELSE IF ( yrel .GT. C_0_R ) THEN; amplify =-C_05_R; ENDIF
     ENDIF

  ELSE
     xi = yrel/thick

     IF      ( iflag .EQ. PROFILE_LINEAR     ) THEN; amplify =-xi
     ELSE IF ( iflag .EQ. PROFILE_TANH       ) THEN; amplify = C_05_R*TANH(-C_05_R*xi)
     ELSE IF ( iflag .EQ. PROFILE_ERF        ) THEN; amplify = C_05_R* ERF(-C_05_R*xi)
     ELSE IF ( iflag .EQ. PROFILE_BICKLEY    ) THEN; amplify = C_1_R/(COSH(C_05_R*xi))**C_2_R
     ELSE IF ( iflag .EQ. PROFILE_GAUSSIAN   ) THEN; amplify = EXP(-C_05_R*xi**C_2_R) 
     ELSE IF ( iflag .EQ. PROFILE_LINEAR_ERF ) THEN; amplify = C_05_R* ERF(-C_05_R*xi)
     ELSE IF ( iflag .EQ. PROFILE_EKMAN_U    ) THEN; amplify = C_1_R -EXP(-xi)*COS(xi)
     ELSE IF ( iflag .EQ. PROFILE_EKMAN_V    ) THEN; amplify =       -EXP(-xi)*SIN(xi)
     ELSE IF ( iflag .EQ. PROFILE_EKMAN_U_P  ) THEN; amplify = C_1_R -EXP(-xi)*COS(xi) ! + perturbation:
        cnought = C_PI_R*C_PI_R/C_4_R/C_4_R       ! Maximum initial Perturbation is at y=pi/2*thick 
        zamp    = SQRT(C_2_R)*xi*EXP( -xi*xi/C_8_R/cnought ) / (thick*thick*C_4_R * cnought)**C_1_5_R
        amplify = amplify + zamp                  ! Add Perturbations
     ELSE IF ( iflag .EQ. PROFILE_PARABOLIC  ) THEN; amplify = (C_1_R+C_05_R*xi)*(C_1_R-C_05_R*xi)
     ENDIF

  ENDIF
  
! mean profile
  FLOW_SHEAR_TEMPORAL = mean + delta*amplify
  
! -------------------------------------------------------------------
! special profiles
! -------------------------------------------------------------------
! two linear-varying layers
  IF ( iflag .EQ. PROFILE_LINEAR_ERF ) THEN
  IF ( yrel .LT. C_0_R ) THEN; FLOW_SHEAR_TEMPORAL = FLOW_SHEAR_TEMPORAL + param(1)*yrel
  ELSE;                        FLOW_SHEAR_TEMPORAL = FLOW_SHEAR_TEMPORAL + param(2)*yrel; ENDIF
  ENDIF
  
! cropped linear
  IF ( iflag .EQ. PROFILE_LINEAR_CROP ) THEN
  IF ( yrel .LT. C_0_R ) THEN; FLOW_SHEAR_TEMPORAL = MIN(param(1)*yrel,param(1)*thick)
  ELSE;                        FLOW_SHEAR_TEMPORAL = MAX(param(2)*yrel,param(2)*thick); ENDIF
  ENDIF

! mixed layer
  IF ( iflag .EQ. PROFILE_MIXEDLAYER ) THEN
  IF ( yrel .LT. C_0_R ) THEN; FLOW_SHEAR_TEMPORAL = MIN(param(1)*yrel,param(1)*thick)
  ELSE;                        FLOW_SHEAR_TEMPORAL = MAX(param(2)*yrel,param(2)*thick); ENDIF
  FLOW_SHEAR_TEMPORAL = FLOW_SHEAR_TEMPORAL - C_025_R*param(2)*thick*(C_1_R -SIGN(C_1_R,y-thick))
  ENDIF

  RETURN
END FUNCTION FLOW_SHEAR_TEMPORAL

