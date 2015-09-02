#include "types.h"

MODULE RAND_LOCAL
  IMPLICIT NONE
  SAVE

! ###################################################################
! Basic options
! ###################################################################
  TINTEGER :: flag_type

! random numbers
  TINTEGER :: isymmetric
  TREAL    :: seed

! spectrum
  TINTEGER :: ispectrum
  TREAL    :: spc_param(5)
  TREAL    :: ucov(6)

! pdf
  TINTEGER :: ipdf

END MODULE RAND_LOCAL
