#include "types.h"

MODULE SCAL_LOCAL
  USE DNS_GLOBAL, ONLY : MAX_NSP

  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_FRC_FREC = 32

! ###################################################################
! Basic options
! ###################################################################
  TINTEGER :: flag_s, flag_mixture

! Geometry and scaling of perturbation
  TREAL    :: thick_ini(MAX_NSP), norm_ini_s(MAX_NSP), ycoor_ini(MAX_NSP)
  TREAL    :: norm_ini_radiation
  
! ###################################################################
! Discrete forcing
! ###################################################################
  TINTEGER :: imode_discrete

! sinusoidal  
  TINTEGER :: nx2d, nx3d, nz3d
  TREAL    :: A2D(MAX_FRC_FREC), Phix2d(MAX_FRC_FREC)
  TREAL    :: A3D(MAX_FRC_FREC), Phix3d(MAX_FRC_FREC), Phiz3d(MAX_FRC_FREC)

! Gaussian
  TREAL    :: delta_discrete

! Alberto Radiation
  TREAL :: rad_ini

END MODULE SCAL_LOCAL
