#include "types.h"

MODULE FLOW_LOCAL
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_FRC_FREC   = 32

! ###################################################################
! Basic options
! ###################################################################
  TINTEGER :: flag_u, flag_t, flag_dilatation, flag_mixture

! Geometry and scaling of perturbation
  TREAL    :: thick_ini, norm_ini_u, norm_ini_p, ycoor_ini

! Boundary conditions
! 0  Free-Slip/Free-Slip
! 1  No-Slip/Free-Slip
! 2  Free-Slip/No-Slip
! 3  No-Slip/No-Slip
  TINTEGER :: flag_wall

! ###################################################################
! Discrete forcing
! ###################################################################
  TINTEGER :: ifrcdsc_mode
  TREAL    :: frc_delta
  
  TINTEGER :: nx2d, nx3d, nz3d
  TREAL    :: A2D(MAX_FRC_FREC), Phix2d(MAX_FRC_FREC)
  TREAL    :: A3D(MAX_FRC_FREC), Phix3d(MAX_FRC_FREC), Phiz3d(MAX_FRC_FREC)

END MODULE FLOW_LOCAL
