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
  USE BOUNDARY_BCS,ONLY: BcsScalJmin, BcsScalJmax
  
  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL, DIMENSION(isize_field)   :: u,v,w, s_is, hs_is
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(jmax,*)        :: wrk1d
  TREAL, DIMENSION(imax,kmax,2)   :: wrk2d

  TREAL, DIMENSION(:), POINTER    :: p_bcs

  TARGET hs_is

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip, ip_b, ip_t, ibc, bcs(2,2)
  TREAL diff


! #######################################################################
  nxy = imax*jmax

  bcs = 0
  
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF


! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero (Dirichlet)
  BcsScalJmin%ref(:,:,is) = C_0_R
  BcsScalJmax%ref(:,:,is) = C_0_R

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
        p_bcs => hs_is(ip_b:); BcsScalJmin%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
     IF ( BcsScalJmax%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
        p_bcs => hs_is(ip_t:); BcsScalJmax%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
     ip_b = ip_b + nxy ! bottom BC address
     ip_t = ip_t + nxy ! top BC address
  ENDDO


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
  IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( BcsScalJmax%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs_is, &
          BcsScalJmin%ref(1,1,is),BcsScalJmax%ref(1,1,is), wrk1d,tmp1,wrk3d)
  ENDIF

  IF ( BcsScalJmin%type(is) .NE. DNS_SFC_STATIC .OR. &
       BcsScalJmax%type(is) .NE. DNS_SFC_STATIC ) THEN
     CALL BOUNDARY_SURFACE_J(is,bcs,s_is,hs_is,tmp1,tmp2,tmp3,wrk1d,wrk2d,wrk3d)
  ENDIF


! -----------------------------------------------------------------------
! Impose BC at Jmin/max
! -----------------------------------------------------------------------
  ip = 1
  ip_t = 1 + imax*(jmax-1)
  DO k = 1,kmax
     hs_is(ip_b:ip_b+imax-1) = BcsScalJmin%ref(1:imax,k,is); ip_b = ip + nxy
     hs_is(ip_t:ip_t+imax-1) = BcsScalJmax%ref(1:imax,k,is); ip_t = ip + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1
