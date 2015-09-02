#include "types.h"
#include "dns_const.h"

SUBROUTINE DNS_PROFILES(y, b_ref, wrk1d)

  USE DNS_GLOBAL, ONLY : MAX_PROF
  USE DNS_GLOBAL, ONLY : jmax, scaley
  USE DNS_GLOBAL, ONLY : iprof_i, ycoor_i, thick_i, delta_i, mean_i, prof_i
  USE DNS_GLOBAL, ONLY : body_param

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(jmax),   INTENT(IN)    :: y
  TREAL, DIMENSION(jmax),   INTENT(OUT)   :: b_ref
  TREAL, DIMENSION(jmax,3), INTENT(INOUT) :: wrk1d
  
! -----------------------------------------------------------------------
  TINTEGER j, flag, iprof
  TREAL ycenter, thick, delta, mean, param(MAX_PROF), FLOW_SHEAR_TEMPORAL
 
! #######################################################################
! Reference state
  ycenter = y(1) + scaley*ycoor_i(1)
  iprof   = iprof_i(1)
  thick   = thick_i(1)
  delta   = delta_i(1)
  mean    = mean_i (1)
  param(:)= prof_i (:,1)
  DO j = 1,jmax
     wrk1d(j,1) = FLOW_SHEAR_TEMPORAL(iprof, thick, delta, mean, ycenter, param, y(j))
     wrk1d(j,3) = C_0_R
  ENDDO
  flag = EQNS_BOD_LINEAR
  CALL FI_BUOYANCY(flag, i1,jmax,i1, body_param, wrk1d(1,1), b_ref, wrk1d(1,3))

  RETURN
END SUBROUTINE DNS_PROFILES
