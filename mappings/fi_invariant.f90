#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2007/09/12 - J.P. Mellado
!#              Created
!#
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

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

  TINTEGER idummy   ! To use old wrappers to calculate derivatives
  TREAL    dummy(1)
! ###################################################################
  idummy = 0

  CALL PARTIAL_X(idummy, nx,ny,nz, idummy, dummy, u, result, i0,i0, dummy, wrk2d,wrk3d)
  CALL PARTIAL_Y(idummy, nx,ny,nz, idummy, dummy, v, tmp1,   i0,i0, dummy, wrk2d,wrk3d)
  result =  result + tmp1
  CALL PARTIAL_Z(idummy, nx,ny,nz, idummy, dummy, w, tmp1,   i0,i0, dummy, wrk2d,wrk3d)
  result =-(result + tmp1)

  RETURN
END SUBROUTINE FI_INVARIANT_P

!########################################################################
! Second invariant of the velocity gradient tensor, Q, as defined
! by Chong et al. (1990). Regions with high positive values of Q are 
! vorticity dominated, regions with high negative values of Q are 
! strain dominated.
! If incompressible, P=0 and then Q=(w^2-2s^2)/4, where P is the first
! invariant, P=-div(v). 
!########################################################################
SUBROUTINE FI_INVARIANT_Q(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, result, tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx,ny,nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: tmp1,tmp2,tmp3
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  result = tmp1*tmp2 + tmp2*tmp3 + tmp3*tmp1

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
! UyVx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result - tmp1*tmp2

! UzWx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result - tmp1*tmp2

! WyVz
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result - tmp1*tmp2

  RETURN
END SUBROUTINE FI_INVARIANT_Q

!########################################################################
! Calculate third invariant of the velocity gradient tensor (-determinant)
! The derivatives v,y and v,z are repeated to avoid more arrays
!########################################################################
SUBROUTINE FI_INVARIANT_R(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, result, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx,ny,nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! Term u,x (v,y w,z-w,y v,z)
! ###################################################################
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  result = tmp5 *( tmp1 *tmp3 -tmp2 *tmp4 )

! ###################################################################
! Term v,x (u,z w,y-w,z u,y)
! ###################################################################
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp5 *( tmp2 *tmp3 -tmp1 *tmp4 )

! ###################################################################
! Term v,x (u,z w,y-w,z u,y)
! ###################################################################
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp5 *( tmp1 *tmp4 -tmp2 *tmp3 )

  result =-result ! set the right sign

  RETURN
END SUBROUTINE FI_INVARIANT_R
