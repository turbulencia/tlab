!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/16 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!# 
!# Coefficients of pentadiagonal 6th-order scheme for 1st derivative
!# 
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE FDM_C1N6M_COEFF()

  USE TLAB_VARS, ONLY : C1N6M_ALPHA, C1N6M_BETA
  USE TLAB_VARS, ONLY : C1N6M_ALPHA2, C1N6M_BETA2
  USE TLAB_VARS, ONLY : C1N6M_A,   C1N6M_B,   C1N6M_C
  USE TLAB_VARS, ONLY : C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
  USE TLAB_VARS, ONLY :            C1N6M_BD2, C1N6M_CD3

  IMPLICIT NONE

! #######################################################################
  ! 10th-order scheme
  ! C1N6M_ALPHA  = c_1_R/C_2_R
  ! C1N6M_BETA   = C_1_R   / C_20_R 
  ! C1N6M_A      = C_17_R  / C_12_R 
  ! C1N6M_B      = C_101_R / C_150_R
  ! C1N6M_C      = C_1_R   / C_100_R
  ! 
  C1N6M_BETA   = (C_2_R/C_5_R )*(-C_1_R/C_3_R +        C1N6M_ALPHA                    )
  C1N6M_A      = (C_1_R/C_6_R )*( C_9_R       +        C1N6M_ALPHA - C_20_R*C1N6M_BETA)
  C1N6M_B      = (C_1_R/C_15_R)*(-C_9_R       + C_32_R*C1N6M_ALPHA + C_62_R*C1N6M_BETA)
  C1N6M_C      = (C_1_R/C_10_R)*( C_1_R       - C_3_R *C1N6M_ALPHA + C_12_R*C1N6M_BETA)
  !  
  C1N6M_AD2    = C1N6M_A / C_2_R 
  C1N6M_BD4    = C1N6M_B / C_4_R
  C1N6M_CD6    = C1N6M_C / C_6_R
  !  
  C1N6M_BD2    = C1N6M_B / C_2_R
  C1N6M_CD3    = C1N6M_C / C_3_R
  !
  C1N6M_ALPHA2 = C_2_R * C1N6M_ALPHA
  C1N6M_BETA2  = C_2_R * C1N6M_BETA

  RETURN
END SUBROUTINE FDM_C1N6M_COEFF