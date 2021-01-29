#include "types.h"
#include "dns_const.h"

FUNCTION FLOW_SHEAR_TEMPORAL(iflag, thick, delta, mean, ycenter, param, y)

  IMPLICIT NONE

  TINTEGER, INTENT(IN) :: iflag                                 ! type of profile
  TREAL,    INTENT(IN) :: thick, delta, mean, ycenter, param(*) ! parameters defining the profile
  TREAL y, FLOW_SHEAR_TEMPORAL

  ! -------------------------------------------------------------------
  TREAL yrel, xi, amplify, zamp, cnought

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
    xi = yrel /thick

    SELECT CASE( iflag )

    CASE( PROFILE_LINEAR )
      amplify =-xi

    CASE( PROFILE_TANH )
      amplify = C_05_R*TANH(-C_05_R*xi)

    CASE( PROFILE_TANH_SYM )
      amplify = C_05_R *( TANH(-C_05_R*(xi-C_05_R*param(5)/thick)) +TANH(C_05_R*(xi+C_05_R*param(5)/thick)) -C_1_R )

    CASE( PROFILE_TANH_ANTISYM )
      amplify = C_025_R*( TANH(-C_05_R*(xi-C_05_R*param(5)/thick)) -TANH(C_05_R*(xi+C_05_R*param(5)/thick)) )

    CASE( PROFILE_ERF,PROFILE_LINEAR_ERF,PROFILE_ERF_ANTISYM,PROFILE_ERF_SURFACE,PROFILE_LINEAR_ERF_SURFACE )
      amplify = C_05_R* ERF(-C_05_R*xi)

    CASE( PROFILE_PARABOLIC )
      amplify = (C_1_R+C_05_R*xi)*(C_1_R-C_05_R*xi)

    CASE( PROFILE_BICKLEY )
      amplify = C_1_R /(COSH(C_05_R*xi))**C_2_R

    CASE( PROFILE_GAUSSIAN )
      amplify = EXP(-C_05_R*xi**C_2_R)

    CASE( PROFILE_GAUSSIAN_SYM )
      amplify = EXP(-C_05_R*(xi-C_05_R*param(5)/thick)**C_2_R) +EXP(C_05_R*(xi+C_05_R*param(5)/thick)**C_2_R)

    CASE( PROFILE_GAUSSIAN_ANTISYM )
      amplify = EXP(-C_05_R*(xi-C_05_R*param(5)/thick)**C_2_R) -EXP(C_05_R*(xi+C_05_R*param(5)/thick)**C_2_R)

    CASE( PROFILE_EKMAN_U )
      amplify = C_1_R -EXP(-xi)*COS(xi)

    CASE( PROFILE_EKMAN_U_P )
      amplify = C_1_R -EXP(-xi)*COS(xi) ! + perturbation:

      cnought = C_PI_R*C_PI_R/C_4_R/C_4_R       ! Maximum initial Perturbation is at y=pi/2*thick
      zamp    = SQRT(C_2_R)*xi*EXP( -xi*xi/C_8_R/cnought ) / (thick*thick*C_4_R * cnought)**C_1_5_R
      amplify = amplify + zamp                  ! Add Perturbations

    CASE( PROFILE_EKMAN_V )
      amplify =       -EXP(-xi)*SIN(xi)

    END SELECT

  ENDIF

  ! mean profile
  IF ( ABS(delta) .GT. C_0_R ) THEN
    FLOW_SHEAR_TEMPORAL = mean + delta *amplify
  ELSE
    FLOW_SHEAR_TEMPORAL = mean
  ENDIF

  ! -------------------------------------------------------------------
  ! special profiles
  ! -------------------------------------------------------------------
  ! two linear-varying layers
  IF ( iflag .EQ. PROFILE_LINEAR_ERF .OR. iflag .EQ. PROFILE_LINEAR_ERF_SURFACE ) THEN
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

  ! adding surface flux
  IF ( iflag .EQ. PROFILE_ERF_SURFACE .OR. iflag .EQ. PROFILE_LINEAR_ERF_SURFACE ) THEN
    xi = y /param(3)
    FLOW_SHEAR_TEMPORAL = FLOW_SHEAR_TEMPORAL + param(4) *C_05_R *( C_1_R +ERF(-C_05_R*xi) )
  ENDIF

  RETURN
END FUNCTION FLOW_SHEAR_TEMPORAL
