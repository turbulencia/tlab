#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Scalar equation, nonlinear term in convective form and the 
!# diffusion term explicit: 3 2nd order + 3 1st order derivatives.
!# BCs need 1 1st order derivatives in Oy
!#
!########################################################################
SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1&
     (is, u,v,w,s_is,hs_is, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : idiffusion, visc, schmidt
  USE BOUNDARY_BCS
  
  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(isize_field)   :: u,v,w, s_is, hs_is
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(jmax,*)        :: wrk1d
  TREAL, DIMENSION(imax,kmax,2)   :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip, ibc, bcs(2,2)
  TREAL diff

! #######################################################################
  nxy = imax*jmax

  bcs = 0
  
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s_is, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s_is, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s_is, tmp4, tmp1, wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,isize_field
     hs_is(ij) = hs_is(ij) + diff*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO
!$omp end do
!$omp end parallel

! #######################################################################
! Boundary conditions
! #######################################################################
  ibc = 0
  wrk2d(:,:,1:2) = C_0_R ! default is dirichlet
  IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( BcsScalJmax%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs_is, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
! Impose bottom at Jmin 
! -----------------------------------------------------------------------
  ip = 1
  DO k = 1,kmax
     hs_is(ip:ip+imax-1) = wrk2d(1:imax,k,1); ip = ip + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top at Jmax
! -----------------------------------------------------------------------
  ip = 1 + imax*(jmax-1) 
  DO k = 1,kmax
     hs_is(ip:ip+imax-1) = wrk2d(1:imax,k,2); ip = ip + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1
