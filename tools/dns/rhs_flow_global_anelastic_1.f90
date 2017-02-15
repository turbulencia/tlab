#include "types.h"
#include "dns_error.h"
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
!# Momentum equations, nonlinear term in convective form and the 
!# viscous term explicit. 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# Scalar needed for the buoyancy term
!#
!########################################################################
!# ARGUMENTS 
!#
!# wrk1d     Out     Reference density profile
!#
!########################################################################
SUBROUTINE  RHS_FLOW_GLOBAL_ANELASTIC_1&
     (dte, u,v,w,h1,h2,h3, z1, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL

IMPLICIT NONE

#include "integers.h"

  TREAL dte
  TREAL, DIMENSION(imax,jmax,kmax) :: u, v, w, h1, h2, h3, z1
  TREAL, DIMENSION(imax,jmax,kmax) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d
  TREAL                            :: wrk1d(isize_wrk1d,*), wrk2d(*)
  TREAL, DIMENSION(jmax)           :: r_ref
  TREAL, DIMENSION(imax,kmax,*)    :: bcs_hb, bcs_ht

! -----------------------------------------------------------------------
  TINTEGER ij, i, j, k

  TREAL dx(1), dy(1), dz(1) ! To use old wrappers to calculate derivatives

! #######################################################################
! Density reference profile
! #######################################################################
  r_ref(:) = C_1_R

! #######################################################################
! Diffusion and convection terms in momentum equations
! #######################################################################
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)

  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
! BCs s.t. derivative d/dy(u) is zero
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp2, i1, i1, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)

  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           h1(i,j,k) = h1(i,j,k) + visc*( tmp6(i,j,k)+tmp5(i,j,k)+tmp4(i,j,k) ) &
                - r_ref(j)*( w(i,j,k)*tmp3(i,j,k)+v(i,j,k)*tmp2(i,j,k)+u(i,j,k)*tmp1(i,j,k) )
        ENDDO
     ENDDO
  ENDDO

  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)

  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
! BCs s.t. to make sure that product vd/dy(v) at the boundary is zero because v is zero. 
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp2, i1, i1, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)

  IF ( buoyancy%type .EQ. EQNS_NONE ) THEN
     DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
        h2(i,j,k) = h2(i,j,k) + visc*( tmp6(i,j,k)+tmp5(i,j,k)+tmp4(i,j,k) ) &
             - r_ref(j)*( w(i,j,k)*tmp3(i,j,k)+v(i,j,k)*tmp2(i,j,k)+u(i,j,k)*tmp1(i,j,k) )
     ENDDO; ENDDO; ENDDO

  ELSE
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, z1, wrk3d, bbackground)
     DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
        h2(i,j,k) = h2(i,j,k) + wrk3d(i,j,k)*buoyancy%vector(2) + visc*( tmp6(i,j,k)+tmp5(i,j,k)+tmp4(i,j,k) ) &
             - r_ref(j)*( w(i,j,k)*tmp3(i,j,k)+v(i,j,k)*tmp2(i,j,k)+u(i,j,k)*tmp1(i,j,k) )
     ENDDO; ENDDO; ENDDO

  ENDIF

  IF ( kmax_total .GT. 1 ) THEN
     CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
          dz, w, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
          dy, w, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
          dx, w, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)

     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
! BCs s.t. derivative d/dy(w) is zero
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, w, tmp2, i1, i1, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, w, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)

     DO k = 1,kmax
        DO j = 1,jmax
           DO i = 1,imax
              h3(i,j,k) = h3(i,j,k) + visc*( tmp6(i,j,k)+tmp5(i,j,k)+tmp4(i,j,k) ) &
                   - r_ref(j)*( w(i,j,k)*tmp3(i,j,k)+v(i,j,k)*tmp2(i,j,k)+u(i,j,k)*tmp1(i,j,k) )
           ENDDO
        ENDDO
     ENDDO

  ENDIF

! #######################################################################
! Pressure term
! #######################################################################
! -----------------------------------------------------------------------
! Poisson equation
! -----------------------------------------------------------------------
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           tmp1(i,j,k) = h1(i,j,k) + r_ref(j)*u(i,j,k)/dte
           tmp2(i,j,k) = h2(i,j,k) + r_ref(j)*v(i,j,k)/dte
           tmp3(i,j,k) = h3(i,j,k) + r_ref(j)*w(i,j,k)/dte
        ENDDO
     ENDDO
  ENDDO
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)

! forcing term in txc2
  DO ij = 1,imax*jmax*kmax
     tmp1(ij,1,1) = tmp6(ij,1,1) + tmp5(ij,1,1) + tmp4(ij,1,1)
  ENDDO

! Neumman BCs s.t. v=0
  DO k = 1,kmax
     DO i = 1,imax
!        tmp1(i,1,   k) = h2(i,1,   k)
!        tmp1(i,jmax,k) = h2(i,jmax,k)
        bcs_hb(i,k,3) = h2(i,1,   k)
        bcs_ht(i,k,3) = h2(i,jmax,k)        
     ENDDO
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
  DO ij = 1,imax*jmax*kmax
     h1(ij,1,1) = h1(ij,1,1) - tmp2(ij,1,1)
     h2(ij,1,1) = h2(ij,1,1) - tmp3(ij,1,1)
     h3(ij,1,1) = h3(ij,1,1) - tmp4(ij,1,1)
  ENDDO

! Impose BCs v=0
  DO k = 1,kmax
     DO i = 1,imax
        h2(i,1,   k) = C_0_R
        h2(i,jmax,k) = C_0_R
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_ANELASTIC_1
