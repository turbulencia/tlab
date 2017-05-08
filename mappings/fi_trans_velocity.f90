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
SUBROUTINE FI_TRANS_VELOCITY(u,v,w, h1,h2,h3, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d)
  
#ifdef USE_OPENMP  
  USE OMP_LIB
#endif
  USE DNS_GLOBAL
  
  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field) :: u, v, w
  TREAL, DIMENSION(isize_field) :: h1, h2, h3
  TREAL, DIMENSION(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d 
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d

! -----------------------------------------------------------------------
  TREAL                         :: dummy_s, dummy_g
  INTEGER                       :: ij, bcs(2,2)

! #######################################################################
  bcs = 0
  
!UX COMPONENT

!Velocity derivatives
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, tmp6, wrk3d, wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, dummy_s, dummy_g )
  dummy_s = stokes
  dummy_g = settling*buoyancy%vector(1)   !Froude number already included in buoyancy%vector
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
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, tmp6, wrk3d, wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, dummy_s, dummy_g )
  dummy_s = stokes
  dummy_g = settling*buoyancy%vector(2)   !Froude number already included in buoyancy%vector
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
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp6, wrk3d, wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, dummy_s, dummy_g )
  dummy_s = stokes
  dummy_g = settling*buoyancy%vector(3)   !Froude number already included in buoyancy%vector
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
