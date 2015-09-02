      SUBROUTINE FILTH2MPPDNU(kmax, ijmax, nx0, nx1, cf, z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Double top-hat filter on a nonuniform mesh
! # Midpoint(trapezoidal) rule
! # Periodic boundary conditions
! # 
! #########################################################

      IMPLICIT NONE

#include "types.h"

      TINTEGER kmax, ijmax
      TINTEGER nx0, nx1
      TREAL cf(nx0+nx1+1,kmax)
      TREAL z1(ijmax,*)
      TREAL zf1(ijmax,*)

      TINTEGER nx

      nx = nx0+nx1

      CALL FILTHMPPDNU(kmax, ijmax, nx, cf, z1, zf1)
      
      RETURN
      END    
