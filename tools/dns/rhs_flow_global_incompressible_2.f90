#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
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
     (dte, u,v,w,h1,h2,h3, s, q,hq, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax
  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type

IMPLICIT NONE

#include "integers.h"

  TREAL dte
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3, s
  TREAL, DIMENSION(isize_field,*) :: q,hq
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d, vaux
  TREAL, DIMENSION(imax,kmax,*)   :: bcs_hb,bcs_ht

  TARGET tmp2, h2

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip_b, ip_t
  TINTEGER ibc
  TREAL alpha

  TREAL, DIMENSION(:), POINTER :: p_bcs

  TREAL dx(1), dy(1), dz(1) ! To use old wrappers to calculate derivatives

! #######################################################################
  nxy    = imax*jmax

! #######################################################################
! Diffusion and convection terms in Ox momentum equations
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp6, i0,i0, i0,i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp5, i0,i0, i0,i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp4, i0,i0, i0,i0, tmp1, wrk1d, wrk2d, wrk3d)

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
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp6, i0,i0, i0,i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp5, i0,i0, i0,i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp4, i0,i0, i0,i0, tmp1, wrk1d, wrk2d, wrk3d)

  DO ij = 1,isize_field
     h2(ij) = h2(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - C_05_R*( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  IF ( kmax_total .GT. 1 ) THEN

  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, w, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, w, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, w, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  
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
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp3, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)
  h1 = h1 - C_05_R*( tmp6 + tmp5 + tmp4 )

  tmp1 = u*v
  tmp2 = v*v
  tmp3 = w*v
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp3, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)
  h2 = h2 - C_05_R*( tmp6 + tmp5 + tmp4 )

  IF ( kmax_total .GT. 1 ) THEN
  tmp1 = u*w
  tmp2 = v*w
  tmp3 = w*w
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp3, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)
  h3 = h3 - C_05_R*( tmp6 + tmp5 + tmp4 )
  ENDIF

! -----------------------------------------------------------------------
! In case the is dilatation
! -----------------------------------------------------------------------
!  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
!       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
!  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
!       dy, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
!       dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
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
  IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(&
          vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
          vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), g(1)%nodes,g(2)%nodes, q,hq)
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
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, tmp1, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc, dy, tmp2, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, tmp3, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)

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
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
       tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i0, dx, tmp1, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, i0, dz, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)

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
