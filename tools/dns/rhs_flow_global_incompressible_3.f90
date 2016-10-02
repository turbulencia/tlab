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
!# Momentum equations, nonlinear term in divergence form and the 
!# viscous term explicit. 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# Scalar needed for the buoyancy term
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE  RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3&
     (dte, x,y,z,dx,dy,dz, u,v,w,h1,h2,h3, s, q,hq, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht,b_ref, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type

IMPLICIT NONE

#include "integers.h"

  TREAL dte
  TREAL, DIMENSION(*)             :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(*)             :: u,v,w, h1,h2,h3, s
  TREAL, DIMENSION(isize_field,*) :: q,hq
  TREAL, DIMENSION(*)             :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d, vaux
  TREAL, DIMENSION(jmax)          :: b_ref
  TREAL, DIMENSION(imax,kmax,*)   :: bcs_hb, bcs_ht

! -----------------------------------------------------------------------
  TINTEGER ij, i, k

! #######################################################################
! Diffusion and convection terms in momentum equations
! #######################################################################
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
  ENDDO

  DO ij = 1,imax*jmax*kmax
     tmp6(ij)=u(ij)*w(ij)
     tmp5(ij)=u(ij)*v(ij)
     tmp4(ij)=u(ij)*u(ij)
  ENDDO
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp6, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp5, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     h1(ij) = h1(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
  ENDDO

! -----------------------------------------------------------------------
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)
  IF ( buoyancy%active(2)  ) THEN
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, b_ref)
     DO ij = 1,imax*jmax*kmax
        h2(ij) = h2(ij) + wrk3d(ij)*buoyancy%vector(2) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
     ENDDO

  ELSE
     DO ij = 1,imax*jmax*kmax
        h2(ij) = h2(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
     ENDDO

  ENDIF

  DO ij = 1,imax*jmax*kmax
     tmp6(ij)=v(ij)*w(ij)
     tmp5(ij)=v(ij)*v(ij)
     tmp4(ij)=v(ij)*u(ij)
  ENDDO
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp6, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
! BCs s.t. to make sure that product vd/dy(v) at the boundary is zero because v is zero. 
!  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
!       dy, tmp5, tmp2, i1, i1, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp5, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  DO ij = 1,imax*jmax*kmax
     h2(ij) = h2(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
  ENDDO

! -----------------------------------------------------------------------
  IF ( kmax_total .GT. 1 ) THEN
     CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
          dz, w, tmp6, i0, i0, i0, i0, tmp3, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
          dy, w, tmp5, i0, i0, i0, i0, tmp2, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
          dx, w, tmp4, i0, i0, i0, i0, tmp1, wrk1d, wrk2d, wrk3d)
     DO ij = 1,imax*jmax*kmax
        h3(ij) = h3(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) )
     ENDDO

     DO ij = 1,imax*jmax*kmax
        tmp6(ij)=w(ij)*w(ij)
        tmp5(ij)=w(ij)*v(ij)
        tmp4(ij)=w(ij)*u(ij)
     ENDDO
     CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
          dz, tmp6, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
          dy, tmp5, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
     CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
          dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
     DO ij = 1,imax*jmax*kmax
        h3(ij) = h3(ij) - ( tmp3(ij) + tmp2(ij) + tmp1(ij) )
     ENDDO
  ENDIF

! -----------------------------------------------------------------------
! Dilatation term for Bcs
! -----------------------------------------------------------------------
!  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
!       dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
!  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
!       dy, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
!  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
!       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
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
  IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(&
          vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
          vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), x,y, q,hq)
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
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp4, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp5, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp6, i0, i0, wrk1d, wrk2d, wrk3d)

! forcing term in txc2
  DO ij = 1,imax*jmax*kmax
     tmp1(ij) = tmp6(ij) + tmp5(ij) + tmp4(ij)
  ENDDO

! Neumman BCs s.t. v=0
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
!        tmp1(ij) = h2(ij)
        bcs_hb(i,k,3) = h2(ij)
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
!        tmp1(ij) = h2(ij)
        bcs_ht(i,k,3) = h2(ij)
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
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, h1, tmp1, i1, i1, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, h3, tmp2, i1, i1, wrk1d, wrk2d, wrk3d)
  DO k = 1,kmax
     DO i = 1,imax
        ij = i                 + imax*jmax*(k-1) ! bottom
        h1(ij) = ( C_4_R*h1(ij+imax) + h1(ij+2*imax) - C_4_R*dy(2)*tmp1(ij+imax) )/C_5_R
        h3(ij) = ( C_4_R*h3(ij+imax) + h3(ij+2*imax) - C_4_R*dy(2)*tmp2(ij+imax) )/C_5_R
        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
        h1(ij) = ( C_4_R*h1(ij-imax) + h1(ij-2*imax) + C_4_R*dy(imax-1)*tmp1(ij-imax) )/C_5_R
        h3(ij) = ( C_4_R*h3(ij-imax) + h3(ij-2*imax) + C_4_R*dy(imax-1)*tmp2(ij-imax) )/C_5_R
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3
