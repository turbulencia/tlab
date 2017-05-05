#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Compute the scalar quantity da/dx_i*du_j/dx_i*da/dx_j (strain1) and
!# that normalized by grad(a)*grad(a) (strain2)
!#
!########################################################################
SUBROUTINE FI_STRAIN_A(nx,ny,nz, a, u,v,w, strain1,strain2, normal1,normal2,normal3, tmp1, wrk2d,wrk3d)

  USE  DNS_GLOBAL, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: a,u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: strain1,strain2
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: normal1,normal2,normal3
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: tmp1 ! Returns grad(a)*grad(a) 
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2), i
  
! ###################################################################
  bcs = 0

! ###################################################################
! Compute normals
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), a, normal1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), a, normal2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), a, normal3, wrk3d, wrk2d,wrk3d)

! Compute gradient terms with u
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = normal1*tmp1*normal1
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal2*tmp1*normal1
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal3*tmp1*normal1
  
! Compute gradient terms with v
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal1*tmp1*normal2
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal2*tmp1*normal2
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal3*tmp1*normal2
  
! Compute gradient terms with 3
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal1*tmp1*normal3
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal2*tmp1*normal3
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1, wrk3d, wrk2d,wrk3d)
  strain1 = strain1 + normal3*tmp1*normal3

! norm of the gradient of a
  tmp1 = normal1**2 + normal2**2 + normal3**2 

! normalize
  DO i = 1,nx*ny*nz
     IF ( tmp1(i) .GT. C_0_R ) THEN
        strain2(i) = strain1(i) /tmp1(i)
     ELSE
        strain2(i) = strain1(i)
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE FI_STRAIN_A
