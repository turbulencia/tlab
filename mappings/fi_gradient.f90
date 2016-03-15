#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2007/09/11 - J.P. Mellado
!#              Created
!# 2016/03/12 - J.P. Mellado
!#              Reducing memory requirement in one array
!#
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
SUBROUTINE FI_GRADIENT(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, s, result, tmp1, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, s,result, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, s,tmp1,   i0,i0, wrk1d,wrk2d,wrk3d)
  result = result*result + tmp1*tmp1
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, s,tmp1,   i0,i0, wrk1d,wrk2d,wrk3d)
  result = result        + tmp1*tmp1

  RETURN
END SUBROUTINE FI_GRADIENT

!########################################################################
! Calculate the scalar gradient production term as given by -(G_i G_j s_ij)
!########################################################################
SUBROUTINE FI_GRADIENT_PRODUCTION(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, s, u,v,w, result, grad_x,grad_y,grad_z, tmp1,tmp2, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s, u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: grad_x,grad_y,grad_z
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! Vorticity vector
! ###################################################################
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, s, grad_x, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, s, grad_y, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, s, grad_z, i0,i0, wrk1d,wrk2d,wrk3d)

! ###################################################################
! Production term
! ###################################################################
! Ux, Vy, Wz
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = tmp1 *grad_x*grad_x + tmp2 *grad_y*grad_y

  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp2 *grad_z*grad_z

! Uy, Vx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result +( tmp1 +tmp2 ) *grad_x*grad_y

! Uz, Wx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result +( tmp1 +tmp2 ) *grad_x*grad_z

! Vz, Wy
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result =-(result + ( tmp1 +tmp2 ) *grad_y*grad_z)

  RETURN
END SUBROUTINE FI_GRADIENT_PRODUCTION

!########################################################################
! Calculate the gradient diffusion term as given by G_i lap G_i
! The diffusivity D is not multiplied here.
!########################################################################
SUBROUTINE FI_GRADIENT_DIFFUSION(iunifx, iunify, iunifz, imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, s, result, grad, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iunifx,iunify,iunifz, imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: grad
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! G_x
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, s, grad, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, grad, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, grad, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, grad, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  result = grad*( tmp1 +tmp2 +tmp3 )

! G_y
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, s, grad, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, grad, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, grad, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, grad, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  result = result + grad*( tmp1+tmp2+tmp3 )

! G_z
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, s, grad, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, grad, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, grad, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, grad, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  result = result + grad*( tmp1 +tmp2 +tmp3 )

  RETURN
END SUBROUTINE FI_GRADIENT_DIFFUSION
