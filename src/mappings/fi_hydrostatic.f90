#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2007/05/08 - J.P. Mellado
!#              Created
!# 2017/03/23 - J.P. Mellado
!#              Rewritten in terms of FDM operators
!#
!########################################################################
!# DESCRIPTION
!#
!# Evaluate the integral \int_pbg%ymean^y dx/H(x), where H(x) is the scale 
!# height in the system
!#
!########################################################################

!########################################################################
! Compute hydrostatic equilibrium from profiles s=(h,q_t).
!########################################################################
SUBROUTINE FI_HYDROSTATIC_H(g, s, e, T,p, wrk1d)

  USE TLAB_TYPES, ONLY : grid_dt

  USE TLAB_VARS, ONLY : imode_eqns
  USE TLAB_VARS, ONLY : pbg, damkohler, buoyancy
  USE THERMO_VARS, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(g%size),   INTENT(IN)    :: e
  TREAL, DIMENSION(g%size),   INTENT(OUT)   :: T,p
  TREAL, DIMENSION(g%size,*), INTENT(INOUT) :: s      ! We calculate equilibrium composition
  TREAL, DIMENSION(g%size,*), INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER iter, niter, ibc, j, jcenter
  TREAL dummy
  
! ###################################################################
! Get the center  
  DO j = 1,g%size
     IF ( g%nodes(j  ) .LE. pbg%ymean .AND. &
          g%nodes(j+1) .GT. pbg%ymean ) THEN
        jcenter = j
        EXIT
     ENDIF
  ENDDO

! Setting the pressure entry to 1 to get 1/RT
  wrk1d(:,6) = C_1_R
  
! Prepare the pentadiagonal system
  ibc = 1                     ! Boundary condition at the bottom for integral calulation 
  CALL INT_C1N6_LHS(g%size, ibc, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
  CALL PENTADFS(g%size-1,        wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))

  niter = 10

  p(:) = pbg%mean             ! initialize iteration
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN       ! Get ql, if necessary
     s(:,3) = C_0_R
  ENDIF
  DO iter = 1,niter           ! iterate
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_DENSITY(i1,g%size,i1, s, e,wrk1d(1,6), wrk1d(1,7))   ! Get 1/RT
        dummy = -C_1_R / SIGN(pbg%parameters(5),buoyancy%vector(2))
     ELSE
        CALL THERMO_AIRWATER_PH_RE(i1,g%size, i1, s(1,2), p, s(1,1), T)
        CALL THERMO_THERMAL_DENSITY(i1,g%size,i1, s(1,2),wrk1d(1,6),T, wrk1d(1,7)) ! Get 1/RT
        dummy = buoyancy%vector(2)
     ENDIF
     wrk1d(:,7) = dummy *wrk1d(:,7)

! Calculate integral
     CALL INT_C1N6_RHS(g%size,i1, ibc, g%jac, wrk1d(1,7), p)
     CALL PENTADSS(g%size-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), p(2))
     p(1) = C_0_R

! Calculate pressure and normalize s.t. p=pbg%mean at y=pbg%ymean_rel
     p(:) = EXP(p(:))
     IF ( ABS(pbg%ymean-g%nodes(jcenter)) .EQ. C_0_R ) THEN
        dummy = p(jcenter)
     ELSE
        dummy = p(jcenter) + (p(jcenter+1)      -p(jcenter)      ) &
                           / (g%nodes(jcenter+1)-g%nodes(jcenter)) *(pbg%ymean-g%nodes(jcenter))
     ENDIF
     dummy = pbg%mean /dummy     
     p(:)  = dummy *p(:)

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Get ql, if necessary
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL THERMO_AIRWATER_PH(i1,g%size,i1, s(1,2),s(1,1), e,p)
        ELSE
           CALL THERMO_AIRWATER_PH_RE(i1,g%size,i1, s(1,2), p, s(1,1), T)
        ENDIF        
     ENDIF
     
  ENDDO

! compute equilibrium values of T
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_TEMPERATURE(i1,g%size,i1, s, e, T)
  ENDIF
  
  RETURN
END SUBROUTINE FI_HYDROSTATIC_H
