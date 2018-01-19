#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculates the residence times for the lagrangian particles
!#
!########################################################################
SUBROUTINE PARTICLE_TIME_RESIDENCE(dtime, particle_number, l_q, l_hq)    

  USE DNS_GLOBAL, ONLY : isize_particle, inb_particle
  USE LAGRANGE_GLOBAL, ONLY: l_y_lambda, l_y_base

  IMPLICIT NONE

  TREAL dtime
  TINTEGER particle_number
  TREAL, DIMENSION(isize_particle,inb_particle) :: l_q, l_hq

  TINTEGER i

  DO i=1,particle_number
     IF (l_q(i,2) .GT. l_y_lambda)THEN
        l_q(i,6)=l_q(i,6)+ dtime   !time cloud droplets spend on cloud-top
     ENDIF
     IF (l_q(i,2) .GT. l_y_base)THEN
        l_hq(i,6)=l_hq(i,6) + dtime   !time cloud droplets spend in intermediate 2/3 of cloud
     ELSEIF (l_q(i,2) .LE. l_y_base)THEN
        l_q(i,6)=C_0_R    !cloud droplets loose memory when "leaving" cloud
        l_hq(i,6)=C_0_R   !cloud droplets loose memory when "leaving" cloud
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE PARTICLE_TIME_RESIDENCE
