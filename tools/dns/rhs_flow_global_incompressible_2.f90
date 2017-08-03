#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Momentum equations, nonlinear term in skew-symmetric form and the 
!# viscous term explicit. 9 2nd order + 9+9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# Scalar needed for the buoyancy term
!#
!########################################################################
SUBROUTINE  RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2&
     (dte, u,v,w,h1,h2,h3, q,hq, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_wrk1d
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax
  USE BOUNDARY_BUFFER

IMPLICIT NONE

  TREAL dte
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,*) :: q,hq
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d
  TREAL, DIMENSION(imax,kmax,*)   :: bcs_hb,bcs_ht

  TARGET tmp2, h2

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip_b, ip_t
  TINTEGER ibc, bcs(2,2)
  TREAL alpha

  TREAL, DIMENSION(:), POINTER :: p_bcs

! #######################################################################
  nxy    = imax*jmax

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! #######################################################################
! Diffusion and convection terms in Ox momentum equations
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), u, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), u, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), u, tmp4, tmp1, wrk2d,wrk3d)

  DO ij = 1,isize_field
     h1(ij) = h1(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - C_05_R*( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO

! BCs s.t. derivative d/dy(u) is zero
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) THEN
     ip_b =                 1
     DO k = 1,kmax
        p_bcs => tmp2(ip_b:); bcs_hb(1:imax,k,1) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     ENDDO
  ENDIF
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN
     ip_t = imax*(jmax-1) + 1
     DO k = 1,kmax
        p_bcs => tmp2(ip_t:); bcs_ht(1:imax,k,1) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
     ENDDO
  ENDIF

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), v, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), v, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), v, tmp4, tmp1, wrk2d,wrk3d)

  DO ij = 1,isize_field
     h2(ij) = h2(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - C_05_R*( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  IF ( g(3)%size .GT. 1 ) THEN

  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), w, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), w, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), w, tmp4, tmp1, wrk2d,wrk3d)
  
  DO ij = 1,isize_field
     h3(ij) = h3(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - C_05_R*( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO

! BCs s.t. derivative d/dy(w) is zero
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) THEN
     ip_b =                 1
     DO k = 1,kmax
        p_bcs => tmp2(ip_b:); bcs_hb(1:imax,k,2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     ENDDO
  ENDIF
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN
     ip_t = imax*(jmax-1) + 1
     DO k = 1,kmax
        p_bcs => tmp2(ip_t:); bcs_ht(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
     ENDDO
  ENDIF

  ENDIF

! #######################################################################
! Adding the conservative-form terms to complete skew-symmetric formulation
! #######################################################################
  tmp1 = u*u
  tmp2 = v*u
  tmp3 = w*u
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d,wrk3d)
  h1 = h1 - C_05_R*( tmp6 + tmp5 + tmp4 )

  tmp1 = u*v
  tmp2 = v*v
  tmp3 = w*v
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d,wrk3d)
  h2 = h2 - C_05_R*( tmp6 + tmp5 + tmp4 )

  IF ( g(3)%size .GT. 1 ) THEN
  tmp1 = u*w
  tmp2 = v*w
  tmp3 = w*w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d,wrk3d)
  h3 = h3 - C_05_R*( tmp6 + tmp5 + tmp4 )
  ENDIF

! -----------------------------------------------------------------------
! In case the is dilatation
! -----------------------------------------------------------------------
!  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
!  DO k = 1,kmax
!     DO i = 1,imax
!        ij = i                 + imax*jmax*(k-1) ! bottom
!        h1(ij) = h1(ij) + C_05_R*u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + C_05_R*w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
!        h1(ij) = h1(ij) + C_05_R*u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + C_05_R*w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     ENDDO
!  ENDDO
!  DO ij = 1,isize_field
!     h1(ij) = h1(ij) - C_05_R*u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     h2(ij) = h2(ij) - C_05_R*v(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     h3(ij) = h3(ij) - C_05_R*w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!  ENDDO 
! no big impact of this part, but in the BCs I need it because I already have
! add 1/2 of it.

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(q, hq)
  ENDIF

! #######################################################################
! Pressure term
! #######################################################################

! -----------------------------------------------------------------------
! Poisson equation
! -----------------------------------------------------------------------
  alpha=C_1_R/dte
  DO ij = 1,isize_field
     tmp1(ij) = h1(ij) + u(ij)*alpha
     tmp2(ij) = h2(ij) + v(ij)*alpha
     tmp3(ij) = h3(ij) + w(ij)*alpha
  ENDDO
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d,wrk3d)

! forcing term in txc2
  DO ij = 1,isize_field
     tmp1(ij) = tmp6(ij) + tmp5(ij) + tmp4(ij)
  ENDDO

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     p_bcs => h2(ip_b:); bcs_hb(1:imax,k,3) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => h2(ip_t:); bcs_ht(1:imax,k,3) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! pressure in tmp1, Oy derivative in tmp3
  ibc = 3
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, ibc, &
       tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient 
! -----------------------------------------------------------------------
  DO ij = 1,isize_field
     h1(ij) = h1(ij) - tmp2(ij)
     h2(ij) = h2(ij) - tmp3(ij)
     h3(ij) = h3(ij) - tmp4(ij)
  ENDDO

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  ibc = 0
  bcs_hb(:,:,1:2) = C_0_R ! default is no-slip (dirichlet)
  bcs_ht(:,:,1:2) = C_0_R
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), h1, bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,tmp1,wrk3d)
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), h3, bcs_hb(1,1,2),bcs_ht(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin 
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     h1(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = C_0_R               ! no penetration
     h3(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2); ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     h1(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = C_0_R               ! no penetration
     h3(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2); ip_t = ip_t + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2
