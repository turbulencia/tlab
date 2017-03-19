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

  TARGET q, tmp3
  
! -----------------------------------------------------------------------
  TINTEGER i, k

  TREAL dx(1), dy(1), dz(1) ! To use old wrappers to calculate derivatives

  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w, rhs
  
! ###################################################################
! Define pointers
  u => q(:,:,:,1)
  v => q(:,:,:,2)
  w => q(:,:,:,3)

! #######################################################################
  p    = C_0_R
  tmp3 = C_0_R
  
  CALL FI_SOURCES_FLOW(u,s, tmp3, wrk1d,wrk3d)
  
! #######################################################################
! Calculate forcing term Ox
! #######################################################################
  rhs => tmp3(:,:,:,1)
  
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, u, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - w *tmp1

  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc,&
       dy, u, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - v *tmp1

  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, u, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - u *tmp1

  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, rhs, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oz
! #######################################################################
  rhs => tmp3(:,:,:,3)
  
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, w, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - w *tmp1

  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc,&
       dy, w, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - v *tmp1

  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, w, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - u *tmp1

  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, rhs, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oy
! #######################################################################
  rhs => tmp3(:,:,:,2)
  
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, v, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - w *tmp1

  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc,&
       dy, v, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - v *tmp1

  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, v, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  rhs = rhs + tmp2 *visc - u *tmp1

! Neumann BCs top and bottom
  DO k = 1,kmax; DO i = 1,imax
     wrk2d(i,k,1) = rhs(i,1   ,k)
     wrk2d(i,k,2) = rhs(i,jmax,k)
  ENDDO; ENDDO
  
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, rhs, tmp1, i0,i0, wrk1d,tmp2,wrk3d)
  p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
  CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, i3, &
       p,wrk3d, tmp1,tmp2, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  RETURN
END SUBROUTINE FI_PRESSURE_BOUSSINESQ
