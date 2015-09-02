! idir_opts(1,*)       number of segments
! idir_opts(2,*) = 1   periodic
!                    0   nonperiodic
! idir_opts(3,*) = 1   mirrored
!                    0   nonmirrored
!
! iseg_opts(1,*,*) = 0   uniform segment
!                    1   Colonius, Lele and Moin stretching
!                    2   Second order polynomial stretching
!                    3   Third order polynomial stretching 
!                    4   Geometric progression
!                    5   Hyperbolic tangent
!                    6   Exponential
#include "types.h"

MODULE GRID_LOCAL
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_PARAMES = 6
  TINTEGER, PARAMETER :: MAX_OPTIONS = 5
  TINTEGER, PARAMETER :: MAX_SEGMENT =10

  TINTEGER, DIMENSION(MAX_OPTIONS,            3) :: idir_opts

  TINTEGER, DIMENSION(            MAX_SEGMENT,3) :: isegdim
  TINTEGER, DIMENSION(MAX_OPTIONS,MAX_SEGMENT,3) :: iseg_opts
  TREAL,    DIMENSION(            MAX_SEGMENT,3) :: isegend
  TREAL,    DIMENSION(MAX_PARAMES,MAX_SEGMENT,3) :: iseg_vals

  TREAL scale(3)

  CHARACTER*36 ifile, ofile, sfile, ffile

END MODULE GRID_LOCAL
