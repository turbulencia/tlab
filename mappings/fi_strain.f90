#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2007/09/11 - J.P. Mellado
!#              Created
!#
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
SUBROUTINE FI_STRAIN_TENSOR(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, result, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),          INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz,6), INTENT(OUT)   :: result ! components Sxx, Syy, Szz, Sxy, Sxz, Syz
  TREAL, DIMENSION(*),          INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! Off diagonal terms; result1 used as auxiliar array
! Uy, Vx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, result(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, result(1,4), i0,i0, wrk1d,wrk2d,wrk3d)
  result(:,4) = C_05_R *( result(:,4) +result(:,1) )

! Uz, Wx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, result(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, result(1,5), i0,i0, wrk1d,wrk2d,wrk3d)
  result(:,5) = C_05_R *( result(:,5) +result(:,1) )

! Vz, Wy
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, result(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, result(1,6), i0,i0, wrk1d,wrk2d,wrk3d)
  result(:,6) = C_05_R *( result(:,6) +result(:,1) )

! diagonal terms
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, result(1,1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, result(1,2), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, result(1,3), i0,i0, wrk1d,wrk2d,wrk3d)

  RETURN
END SUBROUTINE FI_STRAIN_TENSOR

!########################################################################
! Calculate tensor magnitude
!########################################################################
SUBROUTINE FI_STRAIN(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, result, tmp1,tmp2, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! Ux, Vy, Wz
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = tmp1*tmp1 + tmp2*tmp2

  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp2*tmp2

! Uy, Vx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + C_05_R *( tmp1 +tmp2 ) *( tmp1 +tmp2 )

! Uz, Wx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + C_05_R *( tmp1 +tmp2 ) *( tmp1 +tmp2 )

! Vz, Wy
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
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
SUBROUTINE FI_STRAIN_PRODUCTION(imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, result, s_12,s_13,s_23, tmp1,tmp2, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: s_12, s_13, s_23
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! Vorticity production part
! ###################################################################
  CALL FI_VORTICITY_PRODUCTION(nx,ny,nz, u,v,w, result, s_12,s_13,s_23, tmp1,tmp2, wrk2d,wrk3d)
  
  result = C_025_R*result

! ###################################################################
! Pure strain term
! ###################################################################
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  s_12 = C_05_R *( tmp1 +tmp2 )

! Uz, Wx
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  s_13 = C_05_R *( tmp1 +tmp2 )

! Vz, Wy
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  s_23 = C_05_R *( tmp1 +tmp2 )

! -------------------------------------------------------------------
! Production term
! -------------------------------------------------------------------
  result = result + C_2_R *s_12 *s_13 *s_23
     
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp1 *( tmp1*tmp1 +C_3_R *( s_12*s_12 +s_13*s_13 ) )
  
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp1 *( tmp1*tmp1 +C_3_R *( s_12*s_12 +s_23*s_23 ) )
  
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp1 *( tmp1*tmp1 +C_3_R *( s_13*s_13 +s_23*s_23 ) )
  
! right sign
  result =-result
  
  RETURN
END SUBROUTINE FI_STRAIN_PRODUCTION

!########################################################################
!# Calculate the strain diffusion term as given by s_ij lap s_ij
!# The kinematic viscosity \nu is not multiplied here.
!########################################################################
SUBROUTINE FI_STRAIN_DIFFUSION(iunifx,iunify,iunifz, imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w, result, strain, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iunifx, iunify, iunifz, imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: strain
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! -------------------------------------------------------------------
! S_11 = U,x
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, strain, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, strain, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, strain, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, strain, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

  result = strain *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_22 = V,y
! -------------------------------------------------------------------
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, strain, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, strain, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, strain, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, strain, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

  result = result + strain*( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_33 = W,z
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, strain, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, strain, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, strain, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, strain, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

  result = result + strain*( tmp1 +tmp2 +tmp3 )
! -------------------------------------------------------------------
! S_12 = (V,x + U,y)/2
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  strain = tmp1 +tmp2

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, strain, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, strain, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, strain, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

  result = result + C_05_R *strain *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_13 = (U,z + W,x)/2
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  strain = tmp1 +tmp2

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, strain, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, strain, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, strain, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

  result = result + C_05_R *strain *( tmp1 +tmp2 +tmp3 )

! -------------------------------------------------------------------
! S_23 = (W,y + V,z)/2
! -------------------------------------------------------------------
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  strain = tmp1 +tmp2

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, strain, tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, strain, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, strain, tmp1, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

  result = result + C_05_R *strain *( tmp1 +tmp2 +tmp3 )

  RETURN
END SUBROUTINE FI_STRAIN_DIFFUSION

!########################################################################
!# Calculate the strain-pressure term as given by -s_ij p,ij:
!########################################################################
SUBROUTINE FI_STRAIN_PRESSURE(iunifx,iunify,iunifz, imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc, &
     dx,dy,dz, u,v,w,p, result, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iunifx,iunify,iunifz, imode_fdm, nx,ny,nz, i1bc,j1bc,k1bc
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dy,dz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w, p
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! ###################################################################
! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, nx,ny,nz, i1bc, dx, p, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  result = tmp1 *tmp2

  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, nx,ny,nz, j1bc, dy, p, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  result = result +tmp1 *tmp2

  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, nx,ny,nz, k1bc, dz, p, tmp2, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
  result = result +tmp1 *tmp2

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
! 2 s_12 p,xy
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, p, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, tmp2, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, u, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result +tmp1 *( tmp2 +tmp3 )

! 2 s_13 p,xz
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, p, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, tmp2, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, nx,ny,nz, i1bc, dx, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, u, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result +tmp1 *( tmp2 +tmp3 )

! 2 s_23 p,yz
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, p, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, tmp2, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, dy, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, nx,ny,nz, k1bc, dz, v, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
  result = result + tmp1 *( tmp2 +tmp3 )

! right sign
  result =-result

  RETURN
END SUBROUTINE FI_STRAIN_PRESSURE
