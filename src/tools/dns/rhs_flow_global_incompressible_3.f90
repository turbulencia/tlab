#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Momentum equations, nonlinear term in divergence form and the
!# viscous term explicit. 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# Scalar needed for the buoyancy term
!#
!########################################################################
SUBROUTINE  RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3&
     (u,v,w,h1,h2,h3, q,hq, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, isize_wrk1d
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : visc
  USE TIME, only : dte
  USE BOUNDARY_BUFFER
  use OPR_PARTIAL

IMPLICIT NONE

  TREAL, DIMENSION(*)             :: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,*) :: q,hq
  TREAL, DIMENSION(*)             :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(imax,kmax,*)   :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER ij, i, k, ibc, bcs(2,2)

! #######################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! #######################################################################
! Diffusion and convection terms in momentum equations
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), u, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), u, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), u, tmp4, tmp1, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
  ENDDO

  DO ij = 1,imax*jmax*kmax
     tmp6(ij)=u(ij)*w(ij)
     tmp5(ij)=u(ij)*v(ij)
     tmp4(ij)=u(ij)*u(ij)
  ENDDO
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp6, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
  ENDDO

! -----------------------------------------------------------------------
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), v, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), v, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), v, tmp4, tmp1, wrk2d,wrk3d)

  DO ij = 1,imax*jmax*kmax
     h2(ij) = h2(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
  ENDDO

  DO ij = 1,imax*jmax*kmax
     tmp6(ij)=v(ij)*w(ij)
     tmp5(ij)=v(ij)*v(ij)
     tmp4(ij)=v(ij)*u(ij)
  ENDDO
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp6, tmp3, wrk3d, wrk2d,wrk3d)
! BCs s.t. to make sure that product vd/dy(v) at the boundary is zero because v is zero.
!  bcs_loc = 1
!  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_loc, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     h2(ij) = h2(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
  ENDDO

! -----------------------------------------------------------------------
  IF ( g(3)%size .GT. 1 ) THEN
     CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs, g(3), w, tmp6, tmp3, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs, g(2), w, tmp5, tmp2, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs, g(1), w, tmp4, tmp1, wrk2d,wrk3d)
     DO ij = 1,imax*jmax*kmax
        h3(ij) = h3(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
     ENDDO

     DO ij = 1,imax*jmax*kmax
        tmp6(ij)=w(ij)*w(ij)
        tmp5(ij)=w(ij)*v(ij)
        tmp4(ij)=w(ij)*u(ij)
     ENDDO
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp6, tmp3, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
     DO ij = 1,imax*jmax*kmax
        h3(ij) = h3(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
     ENDDO
  ENDIF

! -----------------------------------------------------------------------
! Dilatation term for Bcs
! -----------------------------------------------------------------------
!  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
!  DO k = 1,kmax
!     DO i = 1,imax
!        ij = i                 + imax*jmax*(k-1) ! bottom
!        h1(ij) = h1(ij) + u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
!        h1(ij) = h1(ij) + u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     ENDDO
!  ENDDO

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAX_FLOW()
  ENDIF

! #######################################################################
! Pressure term
! #######################################################################
! -----------------------------------------------------------------------
! Poisson equation
! -----------------------------------------------------------------------
  DO ij = 1,imax*jmax*kmax
     tmp1(ij) = h1(ij) + u(ij)/dte
     tmp2(ij) = h2(ij) + v(ij)/dte
     tmp3(ij) = h3(ij) + w(ij)/dte
  ENDDO
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d,wrk3d)

! forcing term in txc2
  DO ij = 1,imax*jmax*kmax
     tmp1(ij) = tmp6(ij) + tmp5(ij) + tmp4(ij)
  ENDDO

! Neumman BCs s.t. v=0
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
!        tmp1(ij) = h2(ij)
        wrk2d(i,k,1) = h2(ij)
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
!        tmp1(ij) = h2(ij)
        wrk2d(i,k,2) = h2(ij)
     ENDDO
  ENDDO

! pressure in tmp1, Oy derivative in tmp3
  ibc = 3
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, ibc, &
       tmp1,tmp3, tmp2,tmp4, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) - tmp2(ij)
     h2(ij) = h2(ij) - tmp3(ij)
     h3(ij) = h3(ij) - tmp4(ij)
  ENDDO

! #######################################################################
! Boundary conditions
! #######################################################################
! Impose no-penetration BCs v=0
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
        h2(ij) = C_0_R
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
        h2(ij) = C_0_R
     ENDDO
  ENDDO

! Impose free-slip BCs du/dy=0, du/dy=0
  bcs = 1
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), h1, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), h3, tmp2, wrk3d, wrk2d,wrk3d)
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
        h1(ij) = ( C_4_R*h1(ij+imax) + h1(ij+2*imax) - C_4_R*g(2)%jac(2,1)*tmp1(ij+imax) )/C_5_R
        h3(ij) = ( C_4_R*h3(ij+imax) + h3(ij+2*imax) - C_4_R*g(2)%jac(2,1)*tmp2(ij+imax) )/C_5_R
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
        h1(ij) = ( C_4_R*h1(ij-imax) + h1(ij-2*imax) + C_4_R*g(2)%jac(imax-1,1)*tmp1(ij-imax) )/C_5_R
        h3(ij) = ( C_4_R*h3(ij-imax) + h3(ij-2*imax) + C_4_R*g(2)%jac(imax-1,1)*tmp2(ij-imax) )/C_5_R
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3
