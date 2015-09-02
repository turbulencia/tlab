#include "types.h"

MODULE SCAL_LOCAL
  USE DNS_GLOBAL, ONLY : MAX_NSP

  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_DISCRETE = 32

! ###################################################################
! Basic options
! ###################################################################
  TINTEGER :: flag_s, flag_mixture

! Geometry and scaling of perturbation
  TREAL    :: thick_ini(MAX_NSP), norm_ini_s(MAX_NSP), ycoor_ini(MAX_NSP)

! ###################################################################
! Discrete forcing
! ###################################################################
  TINTEGER :: imode_discrete

! sinusoidal  
  TINTEGER :: nx2d, nx3d, nz3d
  TREAL    :: A2D(MAX_DISCRETE), Phix2d(MAX_DISCRETE)
  TREAL    :: A3D(MAX_DISCRETE), Phix3d(MAX_DISCRETE), Phiz3d(MAX_DISCRETE)

! Gaussian
  TREAL    :: delta_discrete

! Alberto Radiation
  TREAL :: rad_ini

END MODULE SCAL_LOCAL
