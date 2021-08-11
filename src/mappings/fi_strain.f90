#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the strain tensor and its magnitude as given by s_ij s_ij:
!# u,x^2 + v,y^2 + w,z^2 +1/2( (u,y+v,x)^2 + (u,z+w,x)^2 + (v,z+w,y)^2 )
!# and terms in its evolution equation.
!#
!########################################################################

!########################################################################
! Calculate tensor components
!########################################################################
SUBROUTINE FI_STRAIN_TENSOR(nx,ny,nz, u,v,w, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6 ! components Sxx, Syy, Szz, Sxy, Sxz, Syz
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! Off diagonal terms; result1 used as auxiliar array
! Uy, Vx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp4, wrk3d, wrk2d,wrk3d)
  tmp4 = C_05_R *( tmp4 +tmp1 )

! Uz, Wx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp5, wrk3d, wrk2d,wrk3d)
  tmp5 = C_05_R *( tmp5 +tmp1 )

! Vz, Wy
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp6, wrk3d, wrk2d,wrk3d)
  tmp6 = C_05_R *( tmp6 +tmp1 )

! diagonal terms
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)

  RETURN
END SUBROUTINE FI_STRAIN_TENSOR

!########################################################################
! Calculate tensor magnitude
!########################################################################
SUBROUTINE FI_STRAIN(nx,ny,nz, u,v,w, result, tmp1,tmp2, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
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

! Ux, Vy, Wz
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = tmp1*tmp1 + tmp2*tmp2

  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + tmp2*tmp2

! Uy, Vx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + C_05_R *( tmp1 +tmp2 ) *( tmp1 +tmp2 )

! Uz, Wx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + C_05_R *( tmp1 +tmp2 ) *( tmp1 +tmp2 )

! Vz, Wy
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = result + C_05_R *( tmp1 +tmp2 ) *( tmp1 +tmp2 )

  RETURN
END SUBROUTINE FI_STRAIN

!########################################################################
!# Calculate the strain production term as given by -s_ij s_jk s_ki 
!# minus 1/4 of the vorticity production. The first term is:
!# s_11^3 + s_22^3 +  s_33^3 +
!# 3 s_11 (s_12^2+s_13^2) + 3 s_22 (s_21^2+s_23^2) + 3 s_33 (s_31^2+s_32^2) +
!# 2 s_12 s_23 s_31
!# or
!# 2 s_12 s_23 s_31 +
!# s_11 ( s_11^2 + 3 ( s_12^2 + s_13^2 ) ) +
!# s_22 ( s_22^2 + 3 ( s_12^2 + s_23^2 ) ) +
!# s_33 ( s_33^2 + 3 ( s_13^2 + s_23^2 ) )
!########################################################################
SUBROUTINE FI_STRAIN_PRODUCTION(nx,ny,nz, u,v,w, result, s_12,s_13,s_23, tmp1,tmp2, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: s_12, s_13, s_23
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! ###################################################################
! Vorticity production part
! ###################################################################
  CALL FI_VORTICITY_PRODUCTION(nx,ny,nz, u,v,w, result, s_12,s_13,s_23, tmp1,tmp2, wrk2d,wrk3d)
  
  result = C_025_R*result

! ###################################################################
! Pure strain term
! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  s_12 = C_05_R *( tmp1 +tmp2 )

! Uz, Wx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)
  s_13 = C_05_R *( tmp1 +tmp2 )

! Vz, Wy
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  s_23 = C_05_R *( tmp1 +tmp2 )

! -------------------------------------------------------------------
! Production term
! -------------------------------------------------------------------
  result = result + C_2_R *s_12 *s_13 *s_23
     
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  result = result + tmp1 *( tmp1*tmp1 +C_3_R *( s_12*s_12 +s_13*s_13 ) )
  
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp1, wrk3d, wrk2d,wrk3d)
  result = result + tmp1 *( tmp1*tmp1 +C_3_R *( s_12*s_12 +s_23*s_23 ) )
  
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1, wrk3d, wrk2d,wrk3d)
  result = result + tmp1 *( tmp1*tmp1 +C_3_R *( s_13*s_13 +s_23*s_23 ) )
  
! right sign
  result =-result
  
  RETURN
END SUBROUTINE FI_STRAIN_PRODUCTION

!########################################################################
!# Calculate the strain diffusion term as given by s_ij lap s_ij
!# The kinematic viscosity \nu is not multiplied here.
!########################################################################
SUBROUTINE FI_STRAIN_DIFFUSION(nx,ny,nz, u,v,w, result, strain, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: strain
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! -------------------------------------------------------------------
! S_11 = U,x
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, strain, wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), strain, tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), strain, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), strain, tmp1, tmp4, wrk2d,wrk3d)

  result = strain *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_22 = V,y
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, strain, wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), strain, tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), strain, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), strain, tmp1, tmp4, wrk2d,wrk3d)

  result = result + strain*( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_33 = W,z
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, strain, wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), strain, tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), strain, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), strain, tmp1, tmp4, wrk2d,wrk3d)

  result = result + strain*( tmp1 +tmp2 +tmp3 )
! -------------------------------------------------------------------
! S_12 = (V,x + U,y)/2
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  strain = tmp1 +tmp2

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), strain, tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), strain, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), strain, tmp1, tmp4, wrk2d,wrk3d)

  result = result + C_05_R *strain *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_13 = (U,z + W,x)/2
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp2, wrk3d, wrk2d,wrk3d)
  strain = tmp1 +tmp2

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), strain, tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), strain, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), strain, tmp1, tmp4, wrk2d,wrk3d)

  result = result + C_05_R *strain *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_23 = (W,y + V,z)/2
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  strain = tmp1 +tmp2

  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), strain, tmp3, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), strain, tmp2, tmp4, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), strain, tmp1, tmp4, wrk2d,wrk3d)

  result = result + C_05_R *strain *( tmp1 +tmp2 +tmp3 )

  RETURN
END SUBROUTINE FI_STRAIN_DIFFUSION

!########################################################################
!# Calculate the strain-pressure term as given by -s_ij p,ij:
!########################################################################
SUBROUTINE FI_STRAIN_PRESSURE(nx,ny,nz, u,v,w,p, result, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w, p
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, nx,ny,nz, bcs, g(1), p, tmp2, tmp4,  wrk2d,wrk3d)
  result = tmp1 *tmp2

  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, nx,ny,nz, bcs, g(2), p, tmp2, tmp4,  wrk2d,wrk3d)
  result = result +tmp1 *tmp2

  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2, nx,ny,nz, bcs, g(3), p, tmp2, tmp4,  wrk2d,wrk3d)
  result = result +tmp1 *tmp2

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
! 2 s_12 p,xy
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), p,    tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), tmp2, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v,    tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u,    tmp3, wrk3d, wrk2d,wrk3d)
  result = result +tmp1 *( tmp2 +tmp3 )

! 2 s_13 p,xz
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), p,    tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w,    tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u,    tmp3, wrk3d, wrk2d,wrk3d)
  result = result +tmp1 *( tmp2 +tmp3 )

! 2 s_23 p,yz
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), p,    tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w,    tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v,    tmp3, wrk3d, wrk2d,wrk3d)
  result = result + tmp1 *( tmp2 +tmp3 )

! right sign
  result =-result

  RETURN
END SUBROUTINE FI_STRAIN_PRESSURE
