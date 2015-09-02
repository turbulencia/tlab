#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/11/12 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# itype      In    Flag indicating type of filter:
!#                  1 Compact 4th order
!#                  2 Explicit 6th order
!#                  3 Explicit 4th order
!#                  4 Explicit ADM
!# tmp        In    Auxilar 3D array of size only used in ADM type
!#
!########################################################################
SUBROUTINE OPR_FILTER_Y(itype, nx,ny,nz, j1bc,bcs_jmin,bcs_jmax, u, cy, tmp, wrk1d,wrk2d,wrk3d)
  
  IMPLICIT NONE

  TINTEGER itype, nx, ny, nz, j1bc, bcs_jmin, bcs_jmax
  TREAL, DIMENSION(*)    :: u, tmp, cy, wrk2d,wrk3d
  TREAL, DIMENSION(ny,*) :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER nxy, nxz

! #######################################################################
  nxy = nx*ny 
  nxz = nx*nz

!Make y-direction the last one
  CALL DNS_TRANSPOSE(u, nxy, nz, nxy, wrk3d, nz)

! Filter
  IF ( itype .EQ. 1 ) THEN
     CALL FILT4C_KERNEL(ny, nxz, wrk3d, u, j1bc, bcs_jmin, bcs_jmax, wrk1d, cy)
     IF ( j1bc .EQ. 0 ) THEN
        CALL TRIDPFS(ny,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
        CALL TRIDPSS(ny,nxz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u, wrk2d)
     ELSE
        CALL TRIDFS(ny,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
        CALL TRIDSS(ny, nxz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), u)
     ENDIF

  ELSE IF ( itype .EQ. 2 ) THEN
     CALL FILT6E_KERNEL(nxz, ny, j1bc, bcs_jmin, bcs_jmax, wrk3d, u)

  ELSE IF ( itype .EQ. 3 ) THEN
     CALL FILT4E_KERNEL(ny, nxz, j1bc, wrk3d, u, cy)

  ELSE IF ( itype .EQ. 4 ) THEN
     CALL FILTADM_KERNEL(ny, nxz, j1bc, wrk3d, u, tmp, cy)

  ENDIF

! traspose back
  CALL DNS_TRANSPOSE(u, nz, nxy, nz, wrk3d, nxy)

  u(1:nxy*nz) = wrk3d(1:nxy*nz)

  RETURN
END SUBROUTINE OPR_FILTER_Y
