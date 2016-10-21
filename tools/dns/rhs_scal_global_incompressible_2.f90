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
!# Scalar equation, nonlinear term in skew-symmetric form and the 
!# diffusion term explicit. 3 2nd order + 3+3 1st order derivatives.
!#
!########################################################################
SUBROUTINE  RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2&
     (is, u,v,w,s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : bcs_scal_jmin, bcs_scal_jmax

  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL, DIMENSION(isize_field) :: u, v, w, s, hs
  TREAL, DIMENSION(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d
  TREAL, DIMENSION(*)           :: wrk1d
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip, ibc
  TREAL diff

  TREAL dx(1), dy(1), dz(1) ! To use old wrappers to calculate derivatives

! #######################################################################
  nxy = imax*jmax

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, s, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, s, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, s, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

  DO ij = 1,isize_field
     hs(ij) = hs(ij) + diff*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - C_05_R*( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO

! #######################################################################
! Adding the conservative-form terms to complete skew-symmetric formulation
! #######################################################################
  tmp1 = u*s
  tmp2 = v*s
  tmp3 = w*s
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp3, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)
  hs = hs - C_05_R*( tmp6 + tmp5 + tmp4 )

! #######################################################################
! Boundary conditions
! #######################################################################
  ibc = 0
  wrk2d(:,:,1:2) = C_0_R ! default is dirichlet
  IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
! Impose bottom at Jmin 
! -----------------------------------------------------------------------
  ip = 1
  DO k = 1,kmax
     hs(ip:ip+imax-1) = wrk2d(1:imax,k,1); ip = ip + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top at Jmax
! -----------------------------------------------------------------------
  ip = 1 + imax*(jmax-1) 
  DO k = 1,kmax
     hs(ip:ip+imax-1) = wrk2d(1:imax,k,2); ip = ip + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2
