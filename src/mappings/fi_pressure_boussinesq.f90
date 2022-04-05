#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the pressure field from a divergence free velocity field
!# and a body force.
!#
!########################################################################
SUBROUTINE FI_PRESSURE_BOUSSINESQ(q,s, p, tmp1,tmp2,tmp, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imax,jmax,kmax
  USE TLAB_VARS, ONLY : isize_wrk1d, isize_field
  USE TLAB_VARS, ONLY : imode_eqns, istagger
  USE TLAB_VARS, ONLY : rbackground

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field,* ), INTENT(IN),    TARGET :: q,s
  TREAL, DIMENSION(isize_field   ), INTENT(OUT)           :: p                
  TREAL, DIMENSION(isize_field   ), INTENT(INOUT)         :: tmp1,tmp2, wrk3d 
  TREAL, DIMENSION(isize_field,3 ), INTENT(INOUT), TARGET :: tmp              
  TREAL, DIMENSION(isize_wrk1d,16), INTENT(INOUT)         :: wrk1d
  TREAL, DIMENSION(imax,kmax,2   ), INTENT(INOUT)         :: wrk2d
! -----------------------------------------------------------------------
  TINTEGER k, bcs(2,2)
  TINTEGER ip_b, ip_t, nxy

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u,v,w
  TREAL, DIMENSION(:),   POINTER :: tmp3,tmp4,tmp5
  TREAL, DIMENSION(:),   POINTER :: p_bcs

! #######################################################################
  nxy  = imax*jmax
  
  bcs  = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  p    = C_0_R
  tmp  = C_0_R

! Define pointers
  u    => q(:,1)
  v    => q(:,2)
  w    => q(:,3)

! #######################################################################
! Sources
  CALL FI_SOURCES_FLOW(q,s, tmp, tmp1, wrk1d,wrk2d,wrk3d)

  tmp3 => tmp(:,1)
  tmp4 => tmp(:,2)
  tmp5 => tmp(:,3)

! Advection and diffusion terms
  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), u,u,u,    p, tmp1, wrk2d,wrk3d) ! store u transposed in tmp1
  tmp3 = tmp3 + p
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, bcs, g(1), v,u,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains u transposed
  tmp4 = tmp4 + p
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, bcs, g(1), w,u,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains u transposed
  tmp5 = tmp5 + p

  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), v,v,v,    p, tmp1, wrk2d,wrk3d) ! store v transposed in tmp1
  tmp4 = tmp4 + p
  CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), u,v,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains v transposed
  tmp3 = tmp3 + p
  CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), w,v,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains v transposed
  tmp5 = tmp5 + p

  CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), w,w,w,    p, tmp1, wrk2d,wrk3d) ! store w transposed in tmp1
  tmp5 = tmp5 + p
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax, bcs, g(3), v,w,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains w transposed
  tmp4 = tmp4 + p
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax, bcs, g(3), u,w,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains w transposed
  tmp3 = tmp3 + p

! Set p-field back to zero
  p = C_0_R

! Calculate forcing term Ox
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
    CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3)
  ENDIF
  IF (istagger .EQ. 1 ) THEN
    CALL OPR_PARTIAL_X(OPR_P1_INT_VP,  imax,jmax,kmax, bcs, g(1), tmp3,tmp2, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Z(OPR_P0_INT_VP,  imax,jmax,kmax, bcs, g(3), tmp2,tmp1, wrk3d, wrk2d,wrk3d)
  ELSE
    CALL OPR_PARTIAL_X(OPR_P1,         imax,jmax,kmax, bcs, g(1), tmp3,tmp1, wrk3d, wrk2d,wrk3d)
  ENDIF
  p = p + tmp1

! Calculate forcing term Oy
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
    CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp4)
  ENDIF
  IF (istagger .EQ. 1 ) THEN
    CALL OPR_PARTIAL_X(OPR_P0_INT_VP,  imax,jmax,kmax, bcs, g(1), tmp4,tmp2, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1,         imax,jmax,kmax, bcs, g(2), tmp2,tmp3, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Z(OPR_P0_INT_VP,  imax,jmax,kmax, bcs, g(3), tmp3,tmp1, wrk3d, wrk2d,wrk3d)
  ELSE
    CALL OPR_PARTIAL_Y(OPR_P1,         imax,jmax,kmax, bcs, g(2), tmp4,tmp1, wrk3d, wrk2d,wrk3d)
  ENDIF
  p = p + tmp1

! Calculate forcing term Oz
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
    CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp5)
  ENDIF
  IF (istagger .EQ. 1 ) THEN
    CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), tmp5,tmp2, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Z(OPR_P1_INT_VP, imax,jmax,kmax, bcs, g(3), tmp2,tmp1, wrk3d, wrk2d,wrk3d)
  ELSE  
    CALL OPR_PARTIAL_Z(OPR_P1,        imax,jmax,kmax, bcs, g(3), tmp5,tmp1, wrk3d, wrk2d,wrk3d)
  ENDIF
  p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration) 
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  IF ( istagger  .EQ. 1 ) THEN ! todo: only need to stagger upper/lower boundary plane, not full h2-array
    CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), tmp4, tmp5, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(3), tmp5, tmp4, wrk3d, wrk2d,wrk3d)
  ENDIF
  DO k = 1,kmax
    p_bcs => tmp4(ip_b:); wrk2d(1:imax,k,1) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
    p_bcs => tmp4(ip_t:); wrk2d(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! Pressure field in p
  CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, i3, &
       p,wrk3d, tmp1,tmp2, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

! Stagger pressure field p back on velocity grid
  IF (istagger .EQ. 1 ) THEN
    CALL OPR_PARTIAL_Z(OPR_P0_INT_PV, imax,jmax,kmax, bcs, g(3), p,    tmp1, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_X(OPR_P0_INT_PV, imax,jmax,kmax, bcs, g(1), tmp1, p,    wrk3d, wrk2d,wrk3d)
  ENDIF

  RETURN
END SUBROUTINE FI_PRESSURE_BOUSSINESQ
