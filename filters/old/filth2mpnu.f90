SUBROUTINE FILTH2MPNU(kmax, ijmax, nx0, nx1, cf, z1, zf1)

! #########################################################
! # FILTER LIBRARY
! #
! # Double top-hat filter on a nonuniform mesh
! # Midpoint(trapezoidal) rule
! # Free boundary conditions
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

  IF ( nx .EQ. 2 ) THEN
     CALL FILTHMPNU2(kmax, ijmax, cf, z1, zf1)
  ELSE IF (nx .EQ. 4 ) THEN
     CALL FILTHMPNU4(kmax, ijmax, cf, z1, zf1)
  ELSE IF (nx .EQ. 6 ) THEN
     CALL FILTHMPNU6(kmax, ijmax, cf, z1, zf1)
  ELSE
     CALL FILTHMPNU(kmax, ijmax, nx, cf, z1, zf1)
  ENDIF

  RETURN
END SUBROUTINE FILTH2MPNU
