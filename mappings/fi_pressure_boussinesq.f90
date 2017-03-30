#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the pressure field from a divergence free velocity field 
!# and a body force.
!# 
!########################################################################
SUBROUTINE FI_PRESSURE_BOUSSINESQ(q,s, p, tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_wrk1d
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_GLOBAL, ONLY : imode_fdm, iunify,iunifx,iunifz, i1bc,j1bc,k1bc

IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q,s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(OUT)   :: p
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: tmp1,tmp2, wrk3d ! larger arrays for the Poisson solver,
  TREAL, DIMENSION(imax,jmax,kmax,3), INTENT(INOUT) :: tmp3             ! but shape (imax,jmax,kmax) is used
  TREAL, DIMENSION(isize_wrk1d,16),   INTENT(INOUT) :: wrk1d            ! to work out forcing term and BCs for p
  TREAL, DIMENSION(imax,kmax,2),      INTENT(INOUT) :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER i, k

  TREAL dx(1), dy(1), dz(1) ! To use old wrappers to calculate derivatives
  
! #######################################################################
  p    = C_0_R
  tmp3 = C_0_R
  
  CALL FI_SOURCES_FLOW(q,s, tmp3, wrk1d,wrk3d)
  
! #######################################################################
! Calculate forcing term Ox
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, q(1,1,1,1), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + tmp2 *visc - q(:,:,:,3) *tmp1

  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc,&
       dy, q(1,1,1,1), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + tmp2 *visc - q(:,:,:,2) *tmp1

  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, q(1,1,1,1), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + tmp2 *visc - q(:,:,:,1) *tmp1

  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp3(1,1,1,1), tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oz
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, q(1,1,1,3), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + tmp2 *visc - q(:,:,:,3) *tmp1

  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc,&
       dy, q(1,1,1,3), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + tmp2 *visc - q(:,:,:,2) *tmp1

  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, q(1,1,1,3), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + tmp2 *visc - q(:,:,:,1) *tmp1

  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp3(1,1,1,3), tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oy
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, q(1,1,1,2), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + tmp2 *visc - q(:,:,:,3) *tmp1

  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc,&
       dy, q(1,1,1,2), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + tmp2 *visc - q(:,:,:,2) *tmp1

  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, q(1,1,1,2), tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + tmp2 *visc - q(:,:,:,1) *tmp1

! Neumann BCs top and bottom
  DO k = 1,kmax; DO i = 1,imax
     wrk2d(i,k,1) = tmp3(i,1   ,k,2)
     wrk2d(i,k,2) = tmp3(i,jmax,k,2)
  ENDDO; ENDDO
  
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp3(1,1,1,2), tmp1, i0,i0, wrk1d,tmp2,wrk3d)
  p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
  CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, i3, &
       p,wrk3d, tmp1,tmp2, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  RETURN
END SUBROUTINE FI_PRESSURE_BOUSSINESQ
