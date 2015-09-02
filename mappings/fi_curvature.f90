!########################################################################
!# Tool/Library FIELDS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/17 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the curvature \kappa of a scalar field z1 and -div(n), where n
!# is the normal vector n = grad(z1)/|grad(z1)|.
!# it is expanded and calculated as 
!# \kappa=-( lap(z1) - n*grad( n*grad(z1) ) )/|grad(z1)| to use the 
!# second-derivative FD operator
!#
!# Array tmp1 contains G\kappa 
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE FI_CURVATURE(iunifx,iunify,iunifz, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
     dx,dy,dz, z1, result, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iunifx, iunify, iunifz
  TINTEGER imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc

  TREAL, DIMENSION(*) :: dx, dy, dz

  TREAL z1(*), result(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*)

  TREAL wrk1d(imax,*)
  TREAL wrk2d(imax,*)
  TREAL wrk3d(imax,jmax,kmax)

! -------------------------------------------------------------------
  TINTEGER ij

! ###################################################################
! -------------------------------------------------------------------
! |grad(z1)|
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, z1, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, z1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, z1, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     tmp4(ij) = SQRT( tmp1(ij)*tmp1(ij) + tmp2(ij)*tmp2(ij) + tmp3(ij)*tmp3(ij) )
  ENDDO

! -------------------------------------------------------------------
! first derivative terms
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, result, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij) = result(ij)*tmp1(ij)
  ENDDO
! tmp1 is now free

  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij) = result(ij) + tmp1(ij)*tmp2(ij)
  ENDDO
! tmp2 is now free

  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp4, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij) = ( result(ij) + tmp2(ij)*tmp3(ij) )/tmp4(ij)
  ENDDO
! tmp3 is now free

! -------------------------------------------------------------------
! second derivative terms and final calculations
! -------------------------------------------------------------------
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, z1, tmp1, i0,i0,i0,i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, z1, tmp2, i0,i0,i0,i0, tmp3, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij) = result(ij) - (tmp1(ij) + tmp2(ij))
  ENDDO
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, z1, tmp1, i0,i0,i0,i0, tmp3, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     tmp1(ij) = result(ij) - tmp1(ij)
     result(ij) = tmp1(ij)/tmp4(ij)
  ENDDO

  RETURN
END SUBROUTINE FI_CURVATURE
