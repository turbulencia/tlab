#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the culr of the vector (u,v,w) in Cartesian coordinates
!#
!########################################################################

SUBROUTINE FI_CURL(nx,ny,nz, u,v,w, wx,wy,wz, tmp, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imode_ibm
  USE IBM_VARS,  ONLY : ibm_partial
  use OPR_PARTIAL
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: wx,wy,wz, tmp
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)

! ###################################################################
  bcs = 0

! IBM   (if .true., OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
  IF ( imode_ibm == 1 ) ibm_partial = .true.

! v,x-u,y
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, wz,  wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp, wrk3d, wrk2d,wrk3d)
  wz = wz - tmp

! u,z-w,x
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, wy,  wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp, wrk3d, wrk2d,wrk3d)
  wy = wy - tmp

! w,y-v,z
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, wx,  wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp, wrk3d, wrk2d,wrk3d)
  wx = wx - tmp

  tmp = wx*wx + wy*wy + wz*wz

  IF ( imode_ibm == 1 ) THEN
    ibm_partial = .false.  
    call IBM_BCS_FIELD(wx)
    call IBM_BCS_FIELD(wy)
    call IBM_BCS_FIELD(wz)
    call IBM_BCS_FIELD(tmp)
  ENDIF

  RETURN
END SUBROUTINE FI_CURL
