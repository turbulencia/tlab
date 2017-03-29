#include "types.h"
#include "dns_const.h"

! Calculate reference profiles for anelastic formulation
SUBROUTINE FI_PROFILES(wrk1d)

  USE DNS_GLOBAL, ONLY : inb_scal, inb_scal_array
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : pbg, sbg, damkohler
  USE DNS_GLOBAL, ONLY : rbackground, bbackground, pbackground, tbackground, epbackground
  USE DNS_GLOBAL, ONLY : buoyancy
  USE THERMO_GLOBAL, ONLY : imixture, GRATIO

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(g(2)%size,*), INTENT(INOUT) :: wrk1d
  
! -----------------------------------------------------------------------
  TINTEGER j, is
  TREAL ycenter, FLOW_SHEAR_TEMPORAL 
 
! #######################################################################
  rbackground = C_1_R ! defaults
  pbackground = C_1_R
  tbackground = C_1_R
  epbackground= C_0_R

! Construct given thermodynamic profiles
  DO is = 1,inb_scal
     ycenter = g(2)%nodes(1) + g(2)%scale *sbg(is)%ymean
     DO j = 1,g(2)%size
        wrk1d(j,is) = FLOW_SHEAR_TEMPORAL(sbg(is)%type, sbg(is)%thick, sbg(is)%delta, sbg(is)%mean, ycenter, sbg(is)%parameters, g(2)%nodes(j))
     ENDDO
  ENDDO
     
  IF ( pbg%parameters(1) .GT. C_0_R ) THEN
! Calculate derived thermodynamic profiles
     epbackground = (g(2)%nodes - g(2)%nodes(1) - g(2)%scale *pbg%ymean) *GRATIO /pbg%parameters(1)
     
     IF ( buoyancy%active(2) ) THEN
!        CALL FI_HYDROSTATIC_H_OLD(g(2)%size, g(2)%nodes, wrk1d, epbackground, tbackground, pbackground, wrk1d(1,4))
        CALL FI_HYDROSTATIC_H(g(2), wrk1d, epbackground, tbackground, pbackground, wrk1d(1,inb_scal_array+1))
     ENDIF
     
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN
        CALL THERMO_AIRWATER_PH(i1,g(2)%size,i1, wrk1d(1,2),wrk1d(1,1), epbackground,pbackground )
        
     ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN 
        CALL THERMO_AIRWATER_LINEAR(i1,g(2)%size,i1, wrk1d, wrk1d(1,inb_scal_array))
        
     ENDIF
     
     CALL THERMO_ANELASTIC_DENSITY(i1,g(2)%size,i1, wrk1d, epbackground,pbackground, rbackground)
     
  ENDIF
  
! Calculate buoyancy profile  
  IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
     CALL THERMO_ANELASTIC_BUOYANCY(i1,g(2)%size,i1, wrk1d, epbackground,pbackground,rbackground, bbackground)
  ELSE
     IF ( buoyancy%active(2) ) THEN ! Check this
        wrk1d(:,inb_scal_array+1) = C_0_R
        CALL FI_BUOYANCY(buoyancy, i1,g(2)%size,i1, wrk1d, bbackground, wrk1d(1,inb_scal_array+1))
     ENDIF
  ENDIF
  
  RETURN
END SUBROUTINE FI_PROFILES