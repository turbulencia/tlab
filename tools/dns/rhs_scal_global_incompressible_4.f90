#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2011/11/28 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Derived from SCAL_FLOW_GLOBAL_INCOMPRESSIBLE_1
!#
!# Scalar equation, nonlinear term in convective form and the 
!# diffusion term explicit: 3 2nd order + 3 1st order derivatives.
!# BCs need 1 1st order derivatives in Oy
!#
!# The array w tmp6 contains the transposed w from 
!# SCAL_FLOW_GLOBAL_INCOMPRESSIBLE_4
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_4&
     (is, dte, dx,dy,dz, u,v,w,s_is,hs_is, s, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

  USE OMP_LIB
  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : bcs_scal_jmin, bcs_scal_jmax

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL dte
  TREAL, DIMENSION(*)             :: dx,dy,dz
  TREAL, DIMENSION(isize_field)   :: u,v,w, s_is, hs_is
  TREAL, DIMENSION(isize_field,*) :: s
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(jmax,*)        :: wrk1d
  TREAL, DIMENSION(imax,kmax,2)   :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip, ibc
  TREAL diff

! #######################################################################
  nxy = imax*jmax

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
  CALL OPR_BURGERS_Z(i1,is, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, s(1,is),w,tmp6, tmp1, i0,i0, i0,i0, tmp2, wrk2d,wrk3d) ! tmp6 contains w transposed
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, s_is, tmp3, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, s_is, tmp5, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Radiation terms
! -----------------------------------------------------------------------
  IF ( is .EQ. irad_scalar ) THEN ! source term in wrk3d
     CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, s(1,is+1), rad_param,&
          wrk1d(1,1), wrk1d(1,2), wrk1d(1,3), wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs_is(ij) = hs_is(ij) + tmp1(ij) + diff*( tmp5(ij) + tmp3(ij) ) &
             - ( v(ij)*tmp2(ij) + u(ij)*tmp4(ij) ) + wrk3d(ij)
     ENDDO
!$omp end do
!$omp end parallel

! -----------------------------------------------------------------------
  ELSE
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs_is(ij) = hs_is(ij) + tmp1(ij) + diff*( tmp5(ij)+tmp3(ij) ) &
             - ( v(ij)*tmp2(ij) + u(ij)*tmp4(ij) )
     ENDDO
!$omp end do
!$omp end parallel

  ENDIF

! #######################################################################
! Boundary conditions
! #######################################################################
  wrk2d(:,:,1:2) = C_0_R ! default is dirichlet

  ibc = 0
  IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, hs_is, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,tmp1,wrk3d)
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
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_4
