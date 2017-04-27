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

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_wrk1d, imode_eqns
  USE DNS_GLOBAL, ONLY : rbackground
  USE DNS_GLOBAL, ONLY : visc

IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q,s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(OUT)   :: p
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: tmp1,tmp2, wrk3d ! larger arrays for the Poisson solver,
  TREAL, DIMENSION(imax,jmax,kmax,3), INTENT(INOUT) :: tmp3             ! but shape (imax,jmax,kmax) is used
  TREAL, DIMENSION(isize_wrk1d,16),   INTENT(INOUT) :: wrk1d            ! to work out forcing term and BCs for p
  TREAL, DIMENSION(imax,kmax,2),      INTENT(INOUT) :: wrk2d

! -----------------------------------------------------------------------
  TINTEGER i, k, bcs(2,2)

! #######################################################################
  bcs  = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  p    = C_0_R
  tmp3 = C_0_R
  
  CALL FI_SOURCES_FLOW(q,s, tmp3, tmp1, wrk1d,wrk2d,wrk3d)
  
! #######################################################################
! Calculate forcing term Ox
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), q(1,1,1,1),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + tmp2 *visc - q(:,:,:,3) *tmp1

  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,1,1,1),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + tmp2 *visc - q(:,:,:,2) *tmp1

  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), q(1,1,1,1),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,1) = tmp3(:,:,:,1) + tmp2 *visc - q(:,:,:,1) *tmp1

  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3(1,1,1,1))
  ENDIF
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp3(1,1,1,1),tmp1, wrk3d, wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oz
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), q(1,1,1,3),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + tmp2 *visc - q(:,:,:,3) *tmp1

  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,1,1,3),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + tmp2 *visc - q(:,:,:,2) *tmp1

  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), q(1,1,1,3),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,3) = tmp3(:,:,:,3) + tmp2 *visc - q(:,:,:,1) *tmp1

  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3(1,1,1,3))
  ENDIF
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3(1,1,1,3),tmp1, wrk3d, wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Calculate forcing term Oy
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), q(1,1,1,2),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + tmp2 *visc - q(:,:,:,3) *tmp1

  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), q(1,1,1,2),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + tmp2 *visc - q(:,:,:,2) *tmp1

  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), q(1,1,1,2),tmp2, tmp1, wrk2d,wrk3d)
  tmp3(:,:,:,2) = tmp3(:,:,:,2) + tmp2 *visc - q(:,:,:,1) *tmp1

! Neumann BCs top and bottom
  DO k = 1,kmax; DO i = 1,imax
     wrk2d(i,k,1) = tmp3(i,1   ,k,2)
     wrk2d(i,k,2) = tmp3(i,jmax,k,2)
  ENDDO; ENDDO

  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     wrk2d(:,:,1) = wrk2d(:,:,1) *rbackground(1)
     wrk2d(:,:,2) = wrk2d(:,:,2) *rbackground(g(2)%size)
     CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3(1,1,1,2))
  ENDIF  
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp3(1,1,1,2),tmp1, wrk3d, wrk2d,wrk3d)
  p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
  CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, i3, &
       p,wrk3d, tmp1,tmp2, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  RETURN
END SUBROUTINE FI_PRESSURE_BOUSSINESQ
