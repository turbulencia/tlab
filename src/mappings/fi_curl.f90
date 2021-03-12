#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the culr of the vector (u,v,w) in Cartesian coordinates
!#
!########################################################################

SUBROUTINE FI_CURL(nx,ny,nz, u,v,w, wx,wy,wz, tmp, wrk2d,wrk3d)

  USE  DNS_GLOBAL, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: wx,wy,wz, tmp
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0
  
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

  RETURN
END SUBROUTINE FI_CURL
