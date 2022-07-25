#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the magnitude of the vorticity as given by w_i w_i:
!# (v,x-u,y)^2 + (u,z-w,x)^2 + (w,y-v,z)^2
!# and terms in its evolution equation.
!#
!########################################################################

!########################################################################
! Calculate the magnitude of the vorticity
!########################################################################
SUBROUTINE FI_VORTICITY(nx,ny,nz, u,v,w, result, tmp1,tmp2, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imode_ibm
  USE IBM_VARS,  ONLY : ibm_partial
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

  ! IBM   (if .true., OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
  IF ( imode_ibm == 1 ) ibm_partial = .true.
  
! v,x-u,y
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = ( tmp1-tmp2 ) *( tmp1-tmp2 )

! u,z-w,x
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp2, wrk3d, wrk2d,wrk3d)
  result = result +( tmp1-tmp2 ) *( tmp1-tmp2 ) 

! w,y-v,z
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = result +( tmp1-tmp2 ) *( tmp1-tmp2 ) 

  IF ( imode_ibm == 1 ) THEN
    ibm_partial = .false.  
    call IBM_BCS_FIELD(result)
  ENDIF

  RETURN
END SUBROUTINE FI_VORTICITY

!########################################################################
! Calculate the vorticity production term as given by w_i w_j s_ij
!########################################################################
SUBROUTINE FI_VORTICITY_PRODUCTION(nx,ny,nz, u,v,w, result, vort_x,vort_y,vort_z, tmp1,tmp2, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: vort_x,vort_y,vort_z
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! ###################################################################
! Vorticity vector
! ###################################################################
! v,x-u,y
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  vort_z = tmp1 -tmp2

! u,z-w,x
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp2, wrk3d, wrk2d,wrk3d)
  vort_y = tmp1 -tmp2

! w,y-v,z
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  vort_x = tmp1 -tmp2
  
! ###################################################################
! Production term
! ###################################################################
! Ux, Vy, Wz
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = tmp1 *vort_x*vort_x +tmp2 *vort_y*vort_y

  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + tmp2 *vort_z*vort_z

! Uy, Vx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + ( tmp1 +tmp2 ) *vort_x *vort_y

! Uz, Wx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + ( tmp1 +tmp2 ) *vort_x *vort_z

! Vz, Wy
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + ( tmp1 +tmp2 ) *vort_y *vort_z

  RETURN
END SUBROUTINE FI_VORTICITY_PRODUCTION

!########################################################################
! Calculate the vorticity diffusion term as given by w_i lap w_i
! The kinematic viscosity \nu is not multiplied here.
!########################################################################
SUBROUTINE FI_VORTICITY_DIFFUSION(nx,ny,nz, u,v,w, result, vort, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: vort
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0
  
! -------------------------------------------------------------------
! W_z = V,x - U,y
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v,   tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u,   tmp2, wrk3d, wrk2d,wrk3d)
  vort = tmp1 -tmp2

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), vort,tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), vort,tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), vort,tmp1, tmp4, wrk2d,wrk3d)
  result = vort *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! W_y = U,z - W,x
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u,   tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w,   tmp2, wrk3d, wrk2d,wrk3d)
  vort = tmp1 -tmp2

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), vort,tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), vort,tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), vort,tmp1, tmp4, wrk2d,wrk3d)
  result = result + vort *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! W_z = W,y - V,z
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w,   tmp1, wrk3d,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v,   tmp2, wrk3d,wrk2d,wrk3d)
  vort = tmp1 -tmp2

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), vort,tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), vort,tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), vort,tmp1, tmp4, wrk2d,wrk3d)
  result = result + vort *( tmp1 +tmp2 +tmp3 )

  RETURN
END SUBROUTINE FI_VORTICITY_DIFFUSION

!########################################################################
! Calculate the baroclinic production of vorticity
! (grad\rho x grad p)/\rho^2
!########################################################################
SUBROUTINE FI_VORTICITY_BAROCLINIC(nx,ny,nz, r,p, result, tmp1,tmp2, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)    :: r,p
  TREAL, DIMENSION(nx*ny*nz,3), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz),   INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),          INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0
  
! ###################################################################
! r,yp,z-r,zp,y
! ###################################################################
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), r, result(1,1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), p, tmp1,        wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), p, result(1,2), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), r, tmp2,        wrk3d, wrk2d,wrk3d)

  result(:,1) = ( result(:,1) *tmp1(:) -result(:,2) *tmp2(:) )/( r(:) *r(:) )

! ###################################################################
! r,zp,x-r,xp,z
! p,z and r,z are in tmp1 and tmp2 respectively
! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), p, result(1,3), wrk3d, wrk2d,wrk3d)
  result(:,2) = result(:,3)*tmp2(:)
  
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), r, tmp2,        wrk3d, wrk2d,wrk3d)
  result(:,2) = ( result(:,2) -tmp1(:) *tmp2(:) )/( r(:) *r(:) )

! ###################################################################
! r,xp,y-r,yp,x
! p,x and r,x are in result3 and tmp2 respectively
! ###################################################################
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), r, tmp1, wrk3d, wrk2d,wrk3d)
  result(:,3) = result(:,3) *tmp1(:)
  
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), p, tmp1, wrk3d, wrk2d,wrk3d)
  result(:,3) = ( tmp1(:) *tmp2(:) -result(:,3) )/( r(:) *r(:) )

  RETURN
END SUBROUTINE FI_VORTICITY_BAROCLINIC
