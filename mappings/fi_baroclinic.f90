!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/06/23 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the baroclinic production of vorticity
!# (grad\rho x grad p)/\rho^2
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE FI__BAROCLINIC(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
     dx, dy, dz, r, p, result, tmp1, tmp2, wrk1d, wrk2d, wrk3d)

  IMPLICIT NONE

#include "types.h"
#include "integers.h"

  TINTEGER imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax)

  TREAL r(*), p(*), result(imax*jmax*kmax,*)
  TREAL tmp1(*), tmp2(*)

  TREAL wrk1d(imax,*)
  TREAL wrk2d(imax,*)
  TREAL wrk3d(imax,jmax,kmax)

! -------------------------------------------------------------------
  TINTEGER ij

! ###################################################################
! r,yp,z-r,zp,y
! ###################################################################
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, r, result(1,1), i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, p, tmp1,        i0, i0, wrk1d, wrk2d, wrk3d)

  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, p, result(ij,2), i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, r, tmp2,         i0, i0, wrk1d, wrk2d, wrk3d)

  DO ij = 1,imax*jmax*kmax
     result(ij,1) = (result(ij,1)*tmp1(ij)-result(ij,2)*tmp2(ij))/(r(ij)*r(ij))
  ENDDO

! ###################################################################
! r,zp,x-r,xp,z
! p,z and r,z are in tmp1 and tmp2 respectively
! ###################################################################
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, p, result(ij,3), i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij,2) = result(ij,3)*tmp2(ij)
  ENDDO

  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, r, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij,2) = (result(ij,2)-tmp1(ij)*tmp2(ij))/(r(ij)*r(ij))
  ENDDO

! ###################################################################
! r,xp,y-r,yp,x
! p,x and r,x are in result3 and tmp2 respectively
! ###################################################################
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, r, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij,3) = result(ij,3)*tmp1(ij)
  ENDDO

  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, p, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     result(ij,3) = (tmp1(ij)*tmp2(ij)-result(ij,3))/(r(ij)*r(ij))
  ENDDO

  RETURN
END SUBROUTINE FI__BAROCLINIC
