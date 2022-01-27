#include "types.h"

MODULE GRID_LOCAL
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: GTYPE_UNIFORM = 0
  TINTEGER, PARAMETER :: GTYPE_TANH    = 5
  TINTEGER, PARAMETER :: GTYPE_EXP     = 6
  
  TINTEGER, PARAMETER :: MAX_PARAMES = 9
  TINTEGER, PARAMETER :: MAX_OPTIONS = 5
  TINTEGER, PARAMETER :: MAX_SEGMENT =10

  TYPE grid_build_dt                          ! Information to construct grid in one direction
    SEQUENCE
    TINTEGER nseg                             ! number of segments in this direction
    LOGICAL mirrored                          ! It true, mirror the grid
    TREAL fixed_scale                         ! If positive, rescale grid to this value
    TINTEGER size(MAX_SEGMENT)                ! number of points in each segment
    TREAL    end(MAX_SEGMENT)                 ! physical end of each segment
    TINTEGER opts(MAX_OPTIONS,MAX_SEGMENT)    ! 0   uniform segment
                                              ! 1   Colonius, Lele and Moin stretching
                                              ! 2   Second order polynomial stretching
                                              ! 3   Third order polynomial stretching
                                              ! 4   Geometric progression
                                              ! 5   Hyperbolic tangent
                                              ! 6   Exponential
    TREAL    vals(MAX_PARAMES,MAX_SEGMENT)
  END TYPE grid_build_dt

  TYPE(grid_build_dt) g_build(3)

END MODULE GRID_LOCAL
