#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the angle between isosurfaces of a and b at each point;
!# the gradient vectors are calculated first and then the angle between
!# them.
!#
!########################################################################
SUBROUTINE FI_ISOSURFACE_ANGLE(nx,ny,nz, a,b, result, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: a,b
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2), ij
  
! ###################################################################
  bcs = 0

  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), a, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), b, tmp2, wrk3d, wrk2d,wrk3d)
  result = tmp1*tmp2
  tmp3 = tmp1*tmp1         
  tmp4 = tmp2*tmp2         

  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), a, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), b, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + tmp1*tmp2
  tmp3 = tmp3 + tmp1*tmp1         
  tmp4 = tmp4 + tmp2*tmp2         
     
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), a, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), b, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + tmp1*tmp2
  tmp3 = tmp3 + tmp1*tmp1         
  tmp4 = tmp4 + tmp2*tmp2

  DO ij = 1,nx*ny*nz
     IF ( tmp3(ij) .GT. C_0_R .AND. tmp4(ij) .GT. C_0_R ) THEN
        result(ij) = result(ij) /SQRT(tmp3(ij)*tmp4(ij))
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE FI_ISOSURFACE_ANGLE
#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the curvature \kappa of a scalar field s and -div(n), where n
!# is the normal vector n = grad(s)/|grad(s)|.
!# it is expanded and calculated as 
!# \kappa=-( lap(s) - n*grad( n*grad(s) ) )/|grad(s)| to use the 
!# second-derivative FD operator
!#
!# Array tmp1 contains G\kappa 
!#
!########################################################################

SUBROUTINE FI_ISOSURFACE_CURVATURE(nx,ny,nz, s, result, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0
  
! -------------------------------------------------------------------
! |grad(s)|
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), s, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), s, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), s, tmp3, wrk3d, wrk2d,wrk3d)
  tmp4 = SQRT( tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3 )
  
! -------------------------------------------------------------------
! first derivative terms
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), tmp4, result, wrk3d, wrk2d,wrk3d)
  result = result *tmp1
! tmp1 is now free
  
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), tmp4, tmp1,   wrk3d, wrk2d,wrk3d)
  result = result + tmp1*tmp2
! tmp2 is now free
  
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), tmp4, tmp2,   wrk3d, wrk2d,wrk3d)
  result = ( result + tmp2*tmp3 )/tmp4
! tmp3 is now free

! -------------------------------------------------------------------
! second derivative terms and final calculations
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), s, tmp1, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), s, tmp2, tmp3, wrk2d,wrk3d)
  result = result - (tmp1 + tmp2)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), s, tmp1, tmp3, wrk2d,wrk3d)
  tmp1 = result - tmp1
  result = tmp1 /tmp4

  RETURN
END SUBROUTINE FI_ISOSURFACE_CURVATURE
