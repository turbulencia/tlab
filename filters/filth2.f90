SUBROUTINE FILTH2(imod, imax,jmax,kmax, kmax_total, &
     iunifx,iunify,iunifz, i1bc,j1bc,k1bc, n0,n1,&
     cf2x,cf2y,cf2z, z1, zf1, wrk3d)

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
#include "types.h"

  IMPLICIT NONE

  TINTEGER imod
  TINTEGER n0, n1
  TINTEGER imax, jmax, kmax, kmax_total
  TINTEGER iunifx, iunify, iunifz
  TINTEGER i1bc, j1bc, k1bc
  TREAL cf2x(*), cf2y(*), cf2z(*)
  TREAL z1(imax,jmax,kmax)
  TREAL zf1(imax,jmax,kmax)
  TREAL wrk3d(imax,jmax,kmax)

  TINTEGER nx0, nx1, ny0, ny1, nz0, nz1

  nx0 = n0
  ny0 = n0
  nz0 = n0
  nx1 = n1
  ny1 = n1
  nz1 = n1

  IF ( kmax_total .EQ. 1 ) THEN
! ###########
! # 2D case #
! ###########

     IF ( imod .EQ. 3 ) THEN
        CALL FILTH2_X(iunifx, i1bc, imax, jmax, kmax, nx0, nx1, cf2x, z1, wrk3d, zf1)
        CALL FILTH2_Y(iunify, j1bc, imax, jmax, kmax, ny0, ny1, cf2y, wrk3d, zf1, wrk3d)
     ELSE IF ( imod .EQ. 1 .OR. imod .EQ. 2 ) THEN
        CALL FILTH2_X(iunifx, i1bc, imax, jmax, kmax, nx0, nx1, cf2x, z1, zf1, wrk3d)
     ENDIF

  ELSE
! ###########
! # 3D case #
! ###########

     IF ( imod .EQ. 2 ) THEN
        CALL FILTH2_Z(iunifz, k1bc, imax, jmax, kmax, nz0, nz1, cf2z, z1, wrk3d, zf1)
     ELSE IF ( imod .EQ. 3 .OR. imod .EQ. 4 ) THEN
        CALL FILTH2_Z(iunifz, k1bc, imax, jmax, kmax, nz0, nz1, cf2z, z1, zf1, wrk3d)
     ENDIF

     IF ( imod .EQ. 3 ) THEN
        CALL FILTH2_X(iunifx, i1bc, imax, jmax, kmax, nx0, nx1, cf2x, zf1, wrk3d, zf1)
        CALL FILTH2_Y(iunify, j1bc, imax, jmax, kmax, ny0, ny1, cf2y, wrk3d, zf1, wrk3d)
     ELSE IF ( imod .EQ. 2 ) THEN
        CALL FILTH2_X(iunifx, i1bc, imax, jmax, kmax, nx0, nx1, cf2x, wrk3d, zf1, wrk3d)
     ELSE IF ( imod .EQ. 1 ) THEN
        CALL FILTH2_X(iunifx, i1bc, imax, jmax, kmax, nx0, nx1, cf2x, z1, zf1, wrk3d)
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE FILTH2
