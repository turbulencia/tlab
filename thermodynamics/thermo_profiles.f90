#include "types.h"
#include "dns_const.h"

! Calculate buoyancy reference state
SUBROUTINE THERMO_PROFILES(wrk1d)

  USE DNS_GLOBAL, ONLY : inb_scal, inb_scal_array
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : iprof_i, ycoor_i, thick_i, delta_i, mean_i, prof_i
  USE DNS_GLOBAL, ONLY : p_init, mean_rho
  USE DNS_GLOBAL, ONLY : rbackground, bbackground, pbackground
  USE DNS_GLOBAL, ONLY : buoyancy
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(g(2)%size,inb_scal_array+1), INTENT(INOUT) :: wrk1d
  
! -----------------------------------------------------------------------
  TINTEGER j, is
  TREAL ycenter, FLOW_SHEAR_TEMPORAL 
 
! #######################################################################
! Density
  rbackground(:) = mean_rho
  
! Buoyancy
  DO is = 1,inb_scal
     ycenter = g(2)%nodes(1) + g(2)%scale*ycoor_i(is)
     DO j = 1,g(2)%size
        wrk1d(j,is) = FLOW_SHEAR_TEMPORAL(iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i(1,is), g(2)%nodes(j))
     ENDDO
  ENDDO
  IF      ( imixture .EQ. MIXT_TYPE_AIRWATER        ) THEN
     bbackground = C_1_R
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
     CALL THERMO_AIRWATER_LINEAR(i1,g(2)%size,i1, wrk1d(1,1), wrk1d(1,inb_scal_array))
  ENDIF
  wrk1d(:,inb_scal_array+1) = C_0_R
  CALL FI_BUOYANCY(buoyancy, i1,g(2)%size,i1, wrk1d(1,1), bbackground, wrk1d(1,inb_scal_array+1))

! pressure
  pbackground(:) = p_init
  
  RETURN
END SUBROUTINE THERMO_PROFILES
