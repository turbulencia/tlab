#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the three invariants of the velocity gradient tensor, div u 
!#
!########################################################################

!########################################################################
! First invariant
!########################################################################
SUBROUTINE FI_INVARIANT_P(nx,ny,nz, u,v,w, result, tmp1, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d
  
! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0
  
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, result, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp1,   wrk3d, wrk2d,wrk3d)
  result =  result + tmp1
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1,   wrk3d, wrk2d,wrk3d)
  result =-(result + tmp1)

  RETURN
END SUBROUTINE FI_INVARIANT_P

!########################################################################
! First invariant on pressure nodes
! (caution: div(u)=0 condition only holds on pressure nodes)
!########################################################################
SUBROUTINE FI_INVARIANT_P_STAG(nx,ny,nz, u,v,w, result, tmp1, tmp2, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1, tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d
  
! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! dudx
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, nx,ny,nz, bcs, g(1), u,    tmp1,   wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, nx,ny,nz, bcs, g(3), tmp1, tmp2,   wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1,        nx,ny,nz, bcs, g(1), tmp2, result, wrk3d, wrk2d,wrk3d)

! dvdy
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, nx,ny,nz, bcs, g(1), v,    tmp1,   wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, nx,ny,nz, bcs, g(3), tmp1, tmp2,   wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1,        nx,ny,nz, bcs, g(2), tmp2, tmp1,   wrk3d, wrk2d,wrk3d)
  
  result =  result + tmp1
  
! dwdz 
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, nx,ny,nz, bcs, g(1), w,    tmp1,   wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, nx,ny,nz, bcs, g(3), tmp1, tmp2,   wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1,        nx,ny,nz, bcs, g(3), tmp2, tmp1,   wrk3d, wrk2d,wrk3d)

  result =-(result + tmp1)

  RETURN
END SUBROUTINE FI_INVARIANT_P_STAG

!########################################################################
! Second invariant of the velocity gradient tensor, Q, as defined
! by Chong et al. (1990). Regions with high positive values of Q are 
! vorticity dominated, regions with high negative values of Q are 
! strain dominated.
! If incompressible, P=0 and then Q=(w^2-2s^2)/4, where P is the first
! invariant, P=-div(v). 
!########################################################################
SUBROUTINE FI_INVARIANT_Q(nx,ny,nz, u,v,w, result, tmp1,tmp2,tmp3, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
  result = tmp1 *tmp2 +tmp2 *tmp3 +tmp3 *tmp1

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
! UyVx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result -tmp1 *tmp2

! UzWx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result -tmp1 *tmp2

! WyVz
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = result -tmp1 *tmp2

  RETURN
END SUBROUTINE FI_INVARIANT_Q

!########################################################################
! Calculate third invariant of the velocity gradient tensor (-determinant)
! The derivatives v,y and v,z are repeated to avoid more arrays
!########################################################################
SUBROUTINE FI_INVARIANT_R(nx,ny,nz, u,v,w, result, tmp1,tmp2,tmp3,tmp4,tmp5, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! ###################################################################
! Term u,x (v,y w,z-w,y v,z)
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp5, wrk3d, wrk2d,wrk3d)
  result = tmp5 *( tmp1 *tmp3 -tmp2 *tmp4 )

! ###################################################################
! Term v,x (u,z w,y-w,z u,y)
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp5, wrk3d, wrk2d,wrk3d)
  result = result + tmp5 *( tmp2 *tmp3 -tmp1 *tmp4 )

! ###################################################################
! Term v,x (u,z w,y-w,z u,y)
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp5, wrk3d, wrk2d,wrk3d)
  result = result + tmp5 *( tmp1 *tmp4 -tmp2 *tmp3 )

  result =-result ! set the right sign

  RETURN
END SUBROUTINE FI_INVARIANT_R
