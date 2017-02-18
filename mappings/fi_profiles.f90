#include "types.h"
#include "dns_const.h"

! Calculate reference profiles for anelastic formulation
SUBROUTINE FI_PROFILES(wrk1d)

  USE DNS_GLOBAL, ONLY : inb_scal, inb_scal_array
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : iprof_i, ycoor_i, thick_i, delta_i, mean_i, prof_i
  USE DNS_GLOBAL, ONLY : p_scale_height, damkohler
  USE DNS_GLOBAL, ONLY : rbackground, bbackground, pbackground, tbackground, epbackground
  USE DNS_GLOBAL, ONLY : buoyancy
  USE THERMO_GLOBAL, ONLY : imixture, GRATIO

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(g(2)%size,inb_scal_array+1), INTENT(INOUT) :: wrk1d
  
! -----------------------------------------------------------------------
  TINTEGER j, is
  TREAL ycenter, FLOW_SHEAR_TEMPORAL 
 
! #######################################################################
  rbackground = C_1_R ! defaults
  pbackground = C_0_R
  tbackground = C_0_R
  epbackground = C_0_R

! Construct given thermodynamic profiles
  DO is = 1,inb_scal
     ycenter = g(2)%nodes(1) + g(2)%scale*ycoor_i(is)
     DO j = 1,g(2)%size
        wrk1d(j,is) = FLOW_SHEAR_TEMPORAL(iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i(1,is), g(2)%nodes(j))
     ENDDO
  ENDDO

! Calculate derived thermodynamic profiles  
  IF      ( imixture .EQ. MIXT_TYPE_AIRWATER        ) THEN
     epbackground = (g(2)%nodes - g(2)%nodes(1) - g(2)%scale*ycoor_i(1)) *GRATIO /p_scale_height

     CALL FI_HYDROSTATIC_AIRWATER_H(g(2)%size, g(2)%nodes, wrk1d(1,2),wrk1d(1,1), epbackground, tbackground, pbackground, wrk1d(1,4))
     IF ( damkohler(3) .LE. C_0_R )  THEN
        CALL THERMO_AIRWATER_PH(i1,g(2)%size,i1, wrk1d(1,2),wrk1d(1,1), epbackground,pbackground )
     ENDIF

     CALL THERMO_AIRWATER_DENSITY(i1,g(2)%size,i1, wrk1d(1,2),wrk1d(1,1), epbackground,pbackground, rbackground)
     
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
     CALL THERMO_AIRWATER_LINEAR(i1,g(2)%size,i1, wrk1d(1,1), wrk1d(1,inb_scal_array))
     
  ENDIF
  
! Calculate buoyancy profile  
  wrk1d(:,inb_scal_array+1) = C_0_R
  CALL FI_BUOYANCY(buoyancy, i1,g(2)%size,i1, wrk1d(1,1), bbackground, wrk1d(1,inb_scal_array+1))
  
  RETURN
END SUBROUTINE FI_PROFILES
