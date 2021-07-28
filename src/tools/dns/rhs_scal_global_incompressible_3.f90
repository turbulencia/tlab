#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Scalar equation, nonlinear term in divergence form and the 
!# diffusion term explicit. 3 2nd order + 3 1st order derivatives.
!#
!########################################################################
SUBROUTINE  RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3&
     (is, u,v,w,s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : idiffusion, visc, schmidt

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(isize_field) :: u,v,w, s, hs
  TREAL, DIMENSION(isize_field) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER ij, i, k, bcs(2,2)
  TREAL diff

! #######################################################################
  bcs = 0

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), s, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), s, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), s, tmp4, tmp1, wrk2d,wrk3d)
  hs = hs + diff*( tmp6+tmp5+tmp4 )

  DO ij = 1,imax*jmax*kmax
     tmp6(ij)=s(ij) *w(ij)
     tmp5(ij)=s(ij) *v(ij)
     tmp4(ij)=s(ij) *u(ij)
  ENDDO
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp6, tmp3, wrk3d, wrk2d,wrk3d)
! which BCs should I use here ?
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
  hs = hs - ( tmp3 + tmp2 + tmp1 )
     
! -----------------------------------------------------------------------
! Dilatation term
! -----------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
!     hs = hs + s*( tmp3 + tmp2 + tmp1 )
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
        hs(ij) = hs(ij) + s(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
        hs(ij) = hs(ij) + s(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_3
