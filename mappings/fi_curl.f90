!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the vorticity as curl v
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE FI_CURL(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, wx,wy,wz, tmp, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc

  TREAL, DIMENSION(*)        :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz) :: u,v,w
  TREAL, DIMENSION(nx*ny*nz) :: wx,wy,wz, tmp
  TREAL, DIMENSION(*)        :: wrk1d,wrk2d,wrk3d

! -------------------------------------------------------------------

! ###################################################################
! v,x-u,y
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, wz,  i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp, i0,i0, wrk1d,wrk2d,wrk3d)
  wz = wz - tmp

! u,z-w,x
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, wy,  i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp, i0,i0, wrk1d,wrk2d,wrk3d)
  wy = wy - tmp

! w,y-v,z
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, wx,  i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp, i0,i0, wrk1d,wrk2d,wrk3d)
  wx = wx - tmp

  tmp = wx*wx + wy*wy + wz*wz

  RETURN
END SUBROUTINE FI_CURL
