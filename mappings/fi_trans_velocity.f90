#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2011/07/13 - A de Lozar
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculates the drift velocity*(-1) from the velocity fiels and the rhs of the equarion
!# u_Drift-(-1) = stokes*(du/dt + u nabla u -g)
!# The three components of the drift velocity are stored in tmp1, tmp2, tmp3
!#
!# Stokes, settling and Freude number must be given
!# It needs 6 temporal arrays which will be used 
!#
!# Question (do it in rhs_flow???)-> we might need 3 extra auxiliary arrays
!#
!########################################################################
SUBROUTINE FI_TRANS_VELOCITY(&
     dx,dy,dz, u,v,w, h1,h2,h3, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)
  
#ifdef USE_OPENMP  
  USE OMP_LIB
#endif
  USE DNS_GLOBAL
  
  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*)           :: dx, dy, dz
  TREAL, DIMENSION(isize_field) :: u, v, w
  TREAL, DIMENSION(isize_field) :: h1, h2, h3
  TREAL, DIMENSION(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d 
  TREAL, DIMENSION(*)           :: wrk1d
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d

! -----------------------------------------------------------------------
  TREAL                         :: dummy_s, dummy_g
  INTEGER                       :: ij

! #######################################################################
!UX COMPONENT

!Velocity derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, u, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, dummy_s, dummy_g )
  dummy_s = stokes
  dummy_g = settling*body_vector(1)   !Froude number already included in body_vector
  IF (dummy_g .NE. C_0_R) THEN
!$omp do
     DO ij = 1,isize_field			
        tmp1(ij) =  dummy_s*( h1(ij) +  u(ij)*tmp4(ij) + v(ij)*tmp5(ij) + w(ij)*tmp6(ij)  - dummy_g )
     ENDDO
!$omp end do
  ELSE
!$omp do
     DO ij = 1,isize_field			
        tmp1(ij) =  dummy_s*( h1(ij) +  u(ij)*tmp4(ij) + v(ij)*tmp5(ij) + w(ij)*tmp6(ij)  )
     ENDDO
!$omp end do
  ENDIF
!$omp end parallel

! #######################################################################
!UY COMPONENT
!Velocity derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, dummy_s, dummy_g )
  dummy_s = stokes
  dummy_g = settling*body_vector(2)   !Froude number already included in body_vector
  IF (dummy_g .NE. C_0_R) THEN
!$omp do
     DO ij = 1,isize_field			
        tmp2(ij) =  dummy_s*( h2(ij) +  u(ij)*tmp4(ij) + v(ij)*tmp5(ij) + w(ij)*tmp6(ij)  - dummy_g )
     ENDDO
!$omp end do
  ELSE
!$omp do
     DO ij = 1,isize_field			
        tmp2(ij) =  dummy_s*( h2(ij) +  u(ij)*tmp4(ij) + v(ij)*tmp5(ij) + w(ij)*tmp6(ij)  )
     ENDDO
!$omp end do
  ENDIF
!$omp end parallel

! #######################################################################
!UZ COMPONENT
!Velocity derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, dummy_s, dummy_g )
  dummy_s = stokes
  dummy_g = settling*body_vector(3)   !Froude number already included in body_vector
  IF (dummy_g .NE. C_0_R) THEN
!$omp do
     DO ij = 1,isize_field			
        tmp3(ij) =  dummy_s*( h3(ij) +  u(ij)*tmp4(ij) + v(ij)*tmp5(ij) + w(ij)*tmp6(ij)  - dummy_g )
     ENDDO
!$omp end do
  ELSE
!$omp do
     DO ij = 1,isize_field			
        tmp3(ij) =  dummy_s*( h3(ij) +  u(ij)*tmp4(ij) + v(ij)*tmp5(ij) + w(ij)*tmp6(ij)  )
     ENDDO
!$omp end do
  ENDIF
!$omp end parallel

  RETURN
END SUBROUTINE FI_TRANS_VELOCITY
