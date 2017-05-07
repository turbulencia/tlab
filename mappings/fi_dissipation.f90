#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculates turbulent dissipation per unit volume \rho \epsilon = \tau'_{ij}u'_{ij}
!# It assumes constant visocsity
!#
!########################################################################
SUBROUTINE FI_DISSIPATION(flag, nx,ny,nz, u,v,w, eps, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : area,visc
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: flag ! 0 for tau_ji  * u_i,j
                                                    ! 1 for tau'_ij * u'_i,j
  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx,ny,nz), INTENT(OUT)   :: eps
  TREAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4, wrk3d
  TREAL, DIMENSION(ny,5),     INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2), j, i1

! ###################################################################
  bcs = 0
  i1 = 1
  
! Diagonal terms
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
  tmp4 = ( tmp1 +tmp2 +tmp3 ) *C_2_R /C_3_R

! 11
  wrk3d = C_2_R *tmp1 -tmp4 ! )*vis
  IF ( flag .EQ. 1 ) THEN
     CALL AVG_IK_V(nx,ny,nz, ny, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,1), wrk1d(1,2), area)
     DO j = 1,ny
        wrk3d(:,j,:)= wrk3d(:,j,:) -wrk1d(j,1)
     ENDDO
  ENDIF
  eps = wrk3d *tmp1
  
! 22
  wrk3d = C_2_R *tmp2 -tmp4 ! )*vis
  IF ( flag .EQ. 1 ) THEN
     CALL AVG_IK_V(nx,ny,nz, ny, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,1), wrk1d(1,2), area)
     CALL AVG_IK_V(nx,ny,nz, ny, v,     g(1)%jac,g(3)%jac, wrk1d(1,3), wrk1d(1,2), area)
     CALL OPR_PARTIAL_Y(OPR_P1, i1,ny,i1, bcs, g(2), wrk1d(1,3),wrk1d(1,2), wrk3d, wrk2d,wrk3d)
     DO j = 1,ny
        wrk3d(:,j,:)= wrk3d(:,j,:) -wrk1d(j,1)
        tmp2(:,j,:) = tmp2(:,j,:)  -wrk1d(j,2)
     ENDDO
  ENDIF
  eps = eps +wrk3d *tmp2

! 33
  wrk3d = C_2_R *tmp3 -tmp4 ! )*vis
  IF ( flag .EQ. 1 ) THEN
     CALL AVG_IK_V(nx,ny,nz, ny, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,1), wrk1d(1,2), area)
     DO j = 1,ny
        wrk3d(:,j,:)= wrk3d(:,j,:) -wrk1d(j,1)
     ENDDO
  ENDIF
  eps = eps +wrk3d *tmp3
     
! Off-diagonal terms
! 12
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp2, wrk3d, wrk2d,wrk3d)

  wrk3d = tmp1 +tmp2 ! )*vis
  IF ( flag .EQ. 1 ) THEN
     CALL AVG_IK_V(nx,ny,nz, ny, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,1), wrk1d(1,2), area)
     CALL AVG_IK_V(nx,ny,nz, ny, u,     g(1)%jac,g(3)%jac, wrk1d(1,3), wrk1d(1,2), area)
     CALL OPR_PARTIAL_Y(OPR_P1, i1,ny,i1, bcs, g(2), wrk1d(1,3),wrk1d(1,2), wrk3d, wrk2d,wrk3d)
     DO j = 1,ny
        wrk3d(:,j,:)= wrk3d(:,j,:) -wrk1d(j,1)
        tmp1(:,j,:) = tmp1(:,j,:)  -wrk1d(j,2)
     ENDDO
  ENDIF
  eps = eps +wrk3d *( tmp1 +tmp2 )

! 13 term
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp2, wrk3d, wrk2d,wrk3d)

  wrk3d = tmp1 +tmp2 ! )*vis
  IF ( flag .EQ. 1 ) THEN
     CALL AVG_IK_V(nx,ny,nz, ny, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,1), wrk1d(1,2), area)
     DO j = 1,ny
        wrk3d(:,j,:)= wrk3d(:,j,:) -wrk1d(j,1)
     ENDDO
  ENDIF
  eps = eps +wrk3d *( tmp1 +tmp2 )

! 23 term
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)

  wrk3d = tmp1 +tmp2 ! )*vis
  IF ( flag .EQ. 1 ) THEN
     CALL AVG_IK_V(nx,ny,nz, ny, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,1), wrk1d(1,2), area)
     CALL AVG_IK_V(nx,ny,nz, ny, w,     g(1)%jac,g(3)%jac, wrk1d(1,3), wrk1d(1,2), area)
     CALL OPR_PARTIAL_Y(OPR_P1, i1,ny,i1, bcs, g(2), wrk1d(1,3),wrk1d(1,2), wrk3d, wrk2d,wrk3d)
     DO j = 1,ny
        wrk3d(:,j,:)= wrk3d(:,j,:) -wrk1d(j,1)
        tmp1(:,j,:) = tmp1(:,j,:)  -wrk1d(j,2)
     ENDDO
  ENDIF
  eps = eps +wrk3d *( tmp1 +tmp2 )

! Final calculation
  eps = eps *visc

  RETURN
END SUBROUTINE FI_DISSIPATION
