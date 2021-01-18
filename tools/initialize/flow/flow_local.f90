#include "types.h"

MODULE FLOW_LOCAL
  USE DNS_TYPES, ONLY : background_dt, discrete_dt
  IMPLICIT NONE
  SAVE

! ###################################################################
! Basic options
! ###################################################################
  TINTEGER :: flag_u, flag_t, flag_dilatation, flag_mixture

  TYPE(background_dt) :: Kini                   ! Geometry of perturbation of initial boundary condition
  TREAL               :: norm_ini_u, norm_ini_p ! Scaling of perturbation
  TYPE(discrete_dt)   :: fp                     ! Discrete perturbation

! Boundary conditions
! 0  Free-Slip/Free-Slip
! 1  No-Slip/Free-Slip
! 2  Free-Slip/No-Slip
! 3  No-Slip/No-Slip
  TINTEGER :: flag_wall

END MODULE FLOW_LOCAL
