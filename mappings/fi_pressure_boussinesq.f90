#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/11/16 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the pressure field from a divergence free velocity field 
!# and a body force.
!# Make sure that the reference profile that you substract to the
!# buoyancy field is the same as in the RHS_FLOW_* routines.
!# 
!########################################################################
!# ARGUMENTS 
!#
!# p, tmp2-3 Aux   Size needed is isize_txc_field for the Poisson
!#                 solver, but the shape (imax,jmax,kmax) is used to
!#                 work out the forcing term and the BCs for p
!# wrk3d     Aux   Size needed is isize_txc_field for the Poisson
!#
!########################################################################
SUBROUTINE FI_PRESSURE_BOUSSINESQ(u,v,w,s, p, tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : MAX_PROF
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, inb_scal, isize_field, isize_wrk1d
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_GLOBAL, ONLY : icoriolis, rotn_vector, rotn_param
  USE DNS_GLOBAL, ONLY : ibodyforce_x,ibodyforce_y,ibodyforce_z, body_param, body_vector
  USE DNS_GLOBAL, ONLY : iprof_i,thick_i,delta_i,mean_i,ycoor_i,prof_i
  USE DNS_GLOBAL, ONLY : imode_fdm, iunify, scaley

IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax), INTENT(IN)    :: u,v,w, s
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: p
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(INOUT) :: tmp1,tmp2,tmp3, wrk3d ! larger arrays
  TREAL, DIMENSION(isize_wrk1d,16), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(imax,kmax,2),    INTENT(INOUT) :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER ij, i, k, flag, iprof
  TREAL ycenter, thick, delta, mean, param(MAX_PROF), FLOW_SHEAR_TEMPORAL
  TREAL dummy, u_geo, w_geo
  TINTEGER iunifx_loc,iunifz_loc, i1bc_loc,j1bc_loc,k1bc_loc

  TREAL, DIMENSION(:), POINTER :: y, dx,dy,dz

! ###################################################################
! Define pointers
                   dx => g(1)%aux(:,1)
  y => g(2)%nodes; dy => g(2)%aux(:,1)
                   dz => g(3)%aux(:,1)

! #######################################################################
  i1bc_loc = 0; iunifx_loc = 0 ! must be periodic and uniform in xOz
  k1bc_loc = 0; iunifz_loc = 0
  j1bc_loc = 1                 ! must be biased in Oy

  p = C_0_R

  u_geo = COS(rotn_param(1))
  w_geo =-SIN(rotn_param(1))

! #######################################################################
! Calculate forcing term Ox
! #######################################################################
! Validation
!  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc_loc, dx, tmp1, p, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i1, iunifz_loc, imode_fdm, imax,jmax,kmax, k1bc_loc,&
       dz, u, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 =        tmp2 *visc - w *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) =                tmp2(ij,1,1)*visc - w(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO
  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc_loc,&
       dy, u, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 = tmp3 + tmp2 *visc - v *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*visc - v(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO
  CALL PARTIAL_XX(i1, iunifx_loc, imode_fdm, imax,jmax,kmax, i1bc_loc,&
       dx, u, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 = tmp3 + tmp2 *visc - u *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*visc - u(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO

  IF ( ibodyforce_x .EQ. EQNS_NONE ) THEN
  ELSE
     wrk1d(:,1) = C_0_R
     CALL FI_BUOYANCY(ibodyforce_x, imax,jmax,kmax, body_param, s, wrk3d, wrk1d)
     tmp3 = tmp3 + body_vector(1) *wrk3d
     ! DO ij = 1,isize_field
     !    tmp3(ij,1,1) = tmp3(ij,1,1) + body_vector(1)*wrk3d(ij,1,1)
     ! ENDDO
     
  ENDIF

  IF ( icoriolis .EQ. EQNS_COR_NORMALIZED ) THEN
     dummy = rotn_vector(2)
     tmp3 = tmp3 + dummy* ( w_geo -w )
     ! DO ij = 1,isize_field
     !    tmp3(ij,1,1) = tmp3(ij,1,1) + dummy*( w_geo-w(ij,1,1) )
     ! ENDDO
  ENDIF

  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc_loc, dx, tmp3, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oz
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz_loc, imode_fdm, imax,jmax,kmax, k1bc_loc,&
       dz, w, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 =        tmp2 *visc - w *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) =                tmp2(ij,1,1)*visc - w(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO
  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc_loc,&
       dy, w, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 = tmp3 + tmp2 *visc - v *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*visc - v(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO
  CALL PARTIAL_XX(i1, iunifx_loc, imode_fdm, imax,jmax,kmax, i1bc_loc,&
       dx, w, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 = tmp3 + tmp2 *visc - u *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*visc - u(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO

  IF ( ibodyforce_z .EQ. EQNS_NONE ) THEN
  ELSE
     wrk1d(:,1) = C_0_R
     CALL FI_BUOYANCY(ibodyforce_z, imax,jmax,kmax, body_param, s, wrk3d, wrk1d)
     tmp3 = tmp3 + body_vector(3) *wrk3d
     ! DO ij = 1,isize_field
     !    tmp3(ij,1,1) = tmp3(ij,1,1) + body_vector(3)*wrk3d(ij,1,1)
     ! ENDDO
     
  ENDIF

  IF ( icoriolis .EQ. EQNS_COR_NORMALIZED ) THEN
     dummy = rotn_vector(2)
     tmp3 = tmp3 + dummy* ( u -u_geo )
     ! DO ij = 1,isize_field
     !    tmp3(ij,1,1) = tmp3(ij,1,1) + dummy*( u(ij,1,1)-u_geo )
     ! ENDDO
  ENDIF

  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc_loc, dz, tmp3, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oy
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz_loc, imode_fdm, imax,jmax,kmax, k1bc_loc,&
       dz, v, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 =        tmp2 *visc - w *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) =                tmp2(ij,1,1)*visc - w(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO
  CALL PARTIAL_YY(i1, iunify,     imode_fdm, imax,jmax,kmax, j1bc_loc,&
       dy, v, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 = tmp3 + tmp2 *visc - v *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*visc - v(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO
  CALL PARTIAL_XX(i1, iunifx_loc, imode_fdm, imax,jmax,kmax, i1bc_loc,&
       dx, v, tmp2, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
  tmp3 = tmp3 + tmp2 *visc - u *tmp1
  ! DO ij = 1,isize_field
  !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*visc - u(ij,1,1)*tmp1(ij,1,1)
  ! ENDDO

! -----------------------------------------------------------------------
! Buoyancy. So far only in the Oy direction. Remember that body_vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( ibodyforce_y .EQ. EQNS_NONE ) THEN
! Validation
!  DO ij = 1,isize_field; tmp3(ij,1,1) = tmp3(ij,1,1) + s(ij,1,1); ENDDO

  ELSE 
! Reference state
     ycenter = y(1) + scaley*ycoor_i(1)
     iprof   = iprof_i(1)
     thick   = thick_i(1)
     delta   = delta_i(1)
     mean    = mean_i (1)
     param(:)= prof_i (:,1)
     DO ij = 1,jmax
        wrk1d(ij,1) = FLOW_SHEAR_TEMPORAL(iprof, thick, delta, mean, ycenter, param, y(ij))
        wrk1d(ij,3) = C_0_R
     ENDDO
     flag = EQNS_BOD_LINEAR
     CALL FI_BUOYANCY(flag,         i1,jmax,i1,     body_param, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))

     CALL FI_BUOYANCY(ibodyforce_y, imax,jmax,kmax, body_param, s,          tmp2,       wrk1d(1,2))
     tmp3 = tmp3 + tmp2 *body_vector(2)
     ! DO ij = 1,isize_field
     !    tmp3(ij,1,1) = tmp3(ij,1,1) + tmp2(ij,1,1)*body_vector(2)
     ! ENDDO

  ENDIF

! Neumann BCs top and bottom
  DO k = 1,kmax; DO i = 1,imax
     wrk2d(i,k,1) = tmp3(i,1   ,k)
     wrk2d(i,k,2) = tmp3(i,jmax,k)
  ENDDO; ENDDO
  
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc_loc, dy, tmp3, tmp1, i0,i0, wrk1d,tmp2,wrk3d)
  p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
  CALL OPR_POISSON_FXZ(imode_fdm,i1,i3, imax,jmax,kmax, &
       y,dx,dy,dz, p,wrk3d, tmp2,tmp3, &
       wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  RETURN
END SUBROUTINE FI_PRESSURE_BOUSSINESQ
