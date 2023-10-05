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
!# Coefficients of pentadiagonal 6th-order scheme for 1st derivative.
!# WARNING: duplicated code, cf. coeff in FDM_C1N6_Jacobian_Penta!
!# 
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE FDM_C1N6M_COEFF()

  USE FDM_PROCS, ONLY : C1N6M_ALPHA, C1N6M_BETA
  USE FDM_PROCS, ONLY : C1N6M_ALPHA2, C1N6M_BETA2
  USE FDM_PROCS, ONLY : C1N6M_A,   C1N6M_B,   C1N6M_C
  USE FDM_PROCS, ONLY : C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
  USE FDM_PROCS, ONLY :            C1N6M_BD2, C1N6M_CD3

  IMPLICIT NONE

! #######################################################################
! -------------------------------------------------------------------
! Pentadiagonal 6th-order scheme (JCP Lele 1992, Eq. 2.1.10) with 
! similar truncation error as 6th-order tridiagonal scheme
! (JCP Lele 1992, Eq. 2.1.7 with alpha=1/3). Value of beta is derived 
! with table 1 (2.1.7 = 2.1.10). Largest stable value with alpha=0.56. 
! (Simulations are unstable for larger alpha values!)
! (values of alpha and beta can be changed, according to Lele...)
! -------------------------------------------------------------------
  C1N6M_ALPHA  = C_56_R/C_100_R
  C1N6M_BETA   = (C_2_R/C_5_R )*(-C_1_R/C_3_R +        C1N6M_ALPHA                    )
  C1N6M_A      = (C_1_R/C_6_R )*( C_9_R       +        C1N6M_ALPHA - C_20_R*C1N6M_BETA)
  C1N6M_B      = (C_1_R/C_15_R)*(-C_9_R       + C_32_R*C1N6M_ALPHA + C_62_R*C1N6M_BETA)
  C1N6M_C      = (C_1_R/C_10_R)*( C_1_R       - C_3_R *C1N6M_ALPHA + C_12_R*C1N6M_BETA)

! #######################################################################
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