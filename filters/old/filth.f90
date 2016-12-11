#include "types.h"

SUBROUTINE FILTH(imod, imax, jmax, kmax, kmax_total, &
     iunifx, iunify, iunifz, i1bc, j1bc, k1bc, n0, &
     cfx, cfy, cfz, z1, zf1, wrk3d)
  
! #######################################################
! # FILTER LIBRARY
! #
! # imod = 1 for 1D filter in X
! # imod = 2 for 2D filter in X and Z
! # imod = 3 for 3D filter
! # imod = 4 for 1D filter in Z
! #
! # Careful with the duplicate entries in the calls !
! #
! #######################################################
  IMPLICIT NONE

  TINTEGER imod
  TINTEGER n0
  TINTEGER imax, jmax, kmax, kmax_total
  TINTEGER iunifx, iunify, iunifz
  TINTEGER i1bc, j1bc, k1bc
  TREAL cfx(*), cfy(*), cfz(*)
  TREAL z1(imax,jmax,kmax)
  TREAL zf1(imax,jmax,kmax)
  TREAL wrk3d(imax,jmax,kmax)

  TINTEGER nx, ny, nz

  nx = n0; ny = n0; nz = n0

  IF ( kmax_total .EQ. 1 ) THEN
! ###########
! # 2D case #
! ###########

     IF ( imod .EQ. 3 ) THEN
        CALL FILTH_X(iunifx, i1bc, imax, jmax, kmax, nx, cfx, z1, wrk3d, zf1)
        CALL FILTH_Y(iunify, j1bc, imax, jmax, kmax, ny, cfy, wrk3d, zf1, wrk3d)
     ELSE IF ( imod .EQ. 1 .OR. imod .EQ. 2 ) THEN
        CALL FILTH_X(iunifx, i1bc, imax, jmax, kmax, nx, cfx, z1, zf1, wrk3d)
     ENDIF

  ELSE
! ###########
! # 3D case #
! ###########

     IF      ( imod .EQ. 2 ) THEN
        CALL FILTH_Z(iunifz, k1bc, imax,jmax,kmax, nz, cfz, z1, wrk3d, zf1)
     ELSE IF ( imod .EQ. 3 .OR. imod .EQ. 4 ) THEN 
        CALL FILTH_Z(iunifz, k1bc, imax,jmax,kmax, nz, cfz, z1, zf1, wrk3d)
     ENDIF

     IF ( imod .EQ. 3 ) THEN
        CALL FILTH_X(iunifx, i1bc, imax,jmax,kmax, nx, cfx, zf1, wrk3d, zf1)
        CALL FILTH_Y(iunify, j1bc, imax,jmax,kmax, ny, cfy, wrk3d, zf1, wrk3d)
     ELSE IF ( imod .EQ. 2 ) THEN
        CALL FILTH_X(iunifx, i1bc, imax,jmax,kmax, nx, cfx, wrk3d, zf1, wrk3d)
     ELSE IF ( imod .EQ. 1 ) THEN
        CALL FILTH_X(iunifx, i1bc, imax,jmax,kmax, nx, cfx, z1, zf1, wrk3d)
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE FILTH
