#include "types.h"
#include "dns_const.h"

!########################################################################
!#
!# Calculates turbulent dissipation per unit volume \rho \epsilon = \tau'_{ij}u'_{ij}
!# It assumes constant visocsity
!#
!########################################################################
SUBROUTINE FI_DISSIPATION(flag, nx,ny,nz, u,v,w, eps, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : area,visc
  USE AVGS, ONLY: AVG_IK_V
    use OPR_PARTIAL
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

! #######################################################################
! Calculate kinetic energy of fluctuating field per unit volume
! #######################################################################
#define rR(j)     wrk1d(j,1)
#define fU(j)     wrk1d(j,2)
#define fV(j)     wrk1d(j,3)
#define fW(j)     wrk1d(j,4)
#define aux(j)    wrk1d(j,5)

SUBROUTINE FI_RTKE(nx,ny,nz, q, wrk1d,wrk3d)

  USE TLAB_VARS, ONLY : imode_eqns
  USE TLAB_VARS, ONLY : g, area, rbackground
  USE AVGS, ONLY: AVG_IK_V

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx,ny,nz,*), INTENT(IN   ) :: q
  TREAL, DIMENSION(ny,5),       INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(nx,ny,nz),   INTENT(INOUT) :: wrk3d

  ! -----------------------------------------------------------------------
  TINTEGER j

  ! #######################################################################
  SELECT CASE ( imode_eqns )
  CASE( DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL)
    CALL AVG_IK_V(nx,ny,nz, ny, q(1,1,1,5), g(1)%jac,g(3)%jac, rR(1), aux(1), area)
    wrk3d = q(:,:,:,5) *q(:,:,:,1)
    CALL AVG_IK_V(nx,ny,nz, ny, wrk3d,      g(1)%jac,g(3)%jac, fU(1), aux(1), area)
    fU(:) = fU(:) /rR(:)
    wrk3d = q(:,:,:,5) *q(:,:,:,2)
    CALL AVG_IK_V(nx,ny,nz, ny, wrk3d,      g(1)%jac,g(3)%jac, fV(1), aux(1), area)
    fV(:) = fV(:) /rR(:)
    wrk3d = q(:,:,:,5) *q(:,:,:,3)
    CALL AVG_IK_V(nx,ny,nz, ny, wrk3d,      g(1)%jac,g(3)%jac, fW(1), aux(1), area)
    fW(:) = fW(:) /rR(:)

  CASE( DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
    rR(:) = rbackground(:)
    CALL AVG_IK_V(nx,ny,nz, ny, q(1,1,1,1), g(1)%jac,g(3)%jac, fU(1), aux(1), area)
    CALL AVG_IK_V(nx,ny,nz, ny, q(1,1,1,2), g(1)%jac,g(3)%jac, fV(1), aux(1), area)
    CALL AVG_IK_V(nx,ny,nz, ny, q(1,1,1,3), g(1)%jac,g(3)%jac, fW(1), aux(1), area)

  END SELECT

  DO j = 1,ny
    wrk3d(:,j,:) = C_05_R *rR(j) *( (q(:,j,:,1)-fU(j))**2 + (q(:,j,:,2)-fV(j))**2 + (q(:,j,:,3)-fW(j))**2 )
  ENDDO

  RETURN
END SUBROUTINE FI_RTKE

!########################################################################
! Reynolds fluctuations of array a
!########################################################################
SUBROUTINE FI_FLUCTUATION_INPLACE(nx,ny,nz, a)
  USE TLAB_VARS, ONLY : g, area
  USE AVGS, ONLY: AVG_IK

  IMPLICIT NONE

  TINTEGER, INTENT(IN)    :: nx,ny,nz
  TREAL,    INTENT(INOUT) :: a(nx,ny,nz)

  ! -------------------------------------------------------------------
  TREAL dummy
  TINTEGER j

  ! ###################################################################
  DO j = 1,ny
    dummy = AVG_IK(nx,ny,nz, j, a, g(1)%jac,g(3)%jac, area)
    a(:,j,:) = a(:,j,:) - dummy
  END DO

  RETURN
END SUBROUTINE FI_FLUCTUATION_INPLACE
