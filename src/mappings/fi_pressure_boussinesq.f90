#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the pressure field from a divergence free velocity field
!# and a body force.
!#
!########################################################################
SUBROUTINE FI_PRESSURE_BOUSSINESQ(q,s, p, tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_wrk1d, imode_eqns
  USE TLAB_VARS, ONLY : rbackground

IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN),   TARGET :: q,s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(OUT)   :: p                ! larger arrays for the Poisson solver,
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: tmp1,tmp2, wrk3d ! but shape (imax,jmax,kmax) is used
  TREAL, DIMENSION(imax,jmax,kmax,3), INTENT(INOUT) :: tmp3             ! to work out forcing term and BCs for p
  TREAL, DIMENSION(isize_wrk1d,16),   INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(imax,kmax,2),      INTENT(INOUT) :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER i, k, bcs(2,2)

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u,v,w

! #######################################################################
  bcs  = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  p    = C_0_R
  tmp3 = C_0_R

! Define pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

! #######################################################################
! Sources
  CALL FI_SOURCES_FLOW(q,s, tmp3, tmp1, wrk1d,wrk2d,wrk3d)

! Advection and diffusion terms
  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), u,u,u,    p, tmp1, wrk2d,wrk3d) ! store u transposed in tmp1
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + p
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, bcs, g(1), v,u,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains u transposed
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + p
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, bcs, g(1), w,u,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains u transposed
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + p

  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), v,v,v,    p, tmp1, wrk2d,wrk3d) ! store v transposed in tmp1
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + p
  CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), u,v,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains v transposed
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + p
  CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), w,v,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains v transposed
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + p

  CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), w,w,w,    p, tmp1, wrk2d,wrk3d) ! store w transposed in tmp1
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + p
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax, bcs, g(3), v,w,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains w transposed
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + p
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax, bcs, g(3), u,w,tmp1, p, tmp2, wrk2d,wrk3d) ! tmp1 contains w transposed
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + p

! Calculate forcing term Ox
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3(1,1,1,1))
  ENDIF
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp3(1,1,1,1),tmp1, wrk3d, wrk2d,wrk3d)
  p = p + tmp1

! Calculate forcing term Oz
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3(1,1,1,3))
  ENDIF
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3(1,1,1,3),tmp1, wrk3d, wrk2d,wrk3d)
  p = p + tmp1

! Calculate forcing term Oy
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3(1,1,1,2))
  ENDIF
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp3(1,1,1,2),tmp1, wrk3d, wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
! Neumann BCs top and bottom
  DO k = 1,kmax; DO i = 1,imax
     wrk2d(i,k,1) = tmp3(i,1   ,k,2)
     wrk2d(i,k,2) = tmp3(i,jmax,k,2)
  ENDDO; ENDDO

  CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, i3, &
       p,wrk3d, tmp1,tmp2, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  RETURN
END SUBROUTINE FI_PRESSURE_BOUSSINESQ
