#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the magnitude of the scalar gradient as given 
!# by G_i G_i, where G_i is ds/dx_i
!# and terms in its evolution equation.
!#
!########################################################################

!########################################################################
! Calculate the magnitude of the scalar gradient
!########################################################################
SUBROUTINE FI_GRADIENT(nx,ny,nz, s, result, tmp1, wrk2d,wrk3d)

  USE  DNS_GLOBAL, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), s,result, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), s,tmp1,   wrk3d ,wrk2d,wrk3d)
  result = result*result + tmp1*tmp1
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), s,tmp1,   wrk3d ,wrk2d,wrk3d)
  result = result        + tmp1*tmp1

  RETURN
END SUBROUTINE FI_GRADIENT

!########################################################################
! Calculate the scalar gradient production term as given by -(G_i G_j s_ij)
!########################################################################
SUBROUTINE FI_GRADIENT_PRODUCTION(nx,ny,nz, s, u,v,w, result, grad_x,grad_y,grad_z, tmp1,tmp2, wrk2d,wrk3d)

  USE  DNS_GLOBAL, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s, u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: grad_x,grad_y,grad_z
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! ###################################################################
! Vorticity vector
! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), s, grad_x, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), s, grad_y, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), s, grad_z, wrk3d ,wrk2d,wrk3d)

! ###################################################################
! Production term
! ###################################################################
! Ux, Vy, Wz
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d ,wrk2d,wrk3d)
  result = tmp1 *grad_x*grad_x + tmp2 *grad_y*grad_y

  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp2, wrk3d ,wrk2d,wrk3d)
  result = result + tmp2 *grad_z*grad_z

! Uy, Vx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d ,wrk2d,wrk3d)
  result = result +( tmp1 +tmp2 ) *grad_x*grad_y

! Uz, Wx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp2, wrk3d ,wrk2d,wrk3d)
  result = result +( tmp1 +tmp2 ) *grad_x*grad_z

! Vz, Wy
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d ,wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d ,wrk2d,wrk3d)
  result =-(result + ( tmp1 +tmp2 ) *grad_y*grad_z)

  RETURN
END SUBROUTINE FI_GRADIENT_PRODUCTION

!########################################################################
! Calculate the gradient diffusion term as given by G_i lap G_i
! The diffusivity D is not multiplied here.
!########################################################################
SUBROUTINE FI_GRADIENT_DIFFUSION(nx,ny,nz, s, result, grad, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE  DNS_GLOBAL, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: grad
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! G_x
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), s,grad, wrk3d ,wrk2d,wrk3d)

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), grad,tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), grad,tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), grad,tmp1, tmp4, wrk2d,wrk3d)
  result = grad*( tmp1 +tmp2 +tmp3 )

! G_y
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), s,grad, wrk3d ,wrk2d,wrk3d)

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), grad,tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), grad,tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), grad,tmp1, tmp4, wrk2d,wrk3d)
  result = result + grad*( tmp1+tmp2+tmp3 )

! G_z
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), s,grad, wrk3d ,wrk2d,wrk3d)

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), grad,tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), grad,tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), grad,tmp1, tmp4, wrk2d,wrk3d)
  result = result + grad*( tmp1 +tmp2 +tmp3 )

  RETURN
END SUBROUTINE FI_GRADIENT_DIFFUSION
