#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Lagrange
!#
!########################################################################
!# HISTORY
!#
!# 2015/03 - L. Muessle
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculates the residence times for the lagrangian particles
!#
!########################################################################
!# ARGUMENTS 
!#
!# 
!#
!########################################################################
SUBROUTINE PARTICLE_TIME_RESIDENCE(dtime, l_q, l_hq)    

  USE DNS_GLOBAL, ONLY : inb_particle, isize_particle
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI, ONLY : particle_vector, ims_pro
#endif

  IMPLICIT NONE

  TREAL dtime
  TREAL, DIMENSION(isize_particle,inb_particle) :: l_q, l_hq
  TINTEGER l_i, local_isize_particle

#ifdef USE_MPI
   local_isize_particle = particle_vector(ims_pro+1)
#else
   local_isize_particle = particle_number
#endif

      IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
         DO l_i=1,local_isize_particle
            IF (l_q(l_i,2) .GT. l_y_lambda)THEN
               l_q(l_i,6)=l_q(l_i,6)+ dtime   !time cloud droplets spend on cloud-top
            ENDIF
            IF (l_q(l_i,2) .GT. l_y_base)THEN
               l_hq(l_i,6)=l_hq(l_i,6) + dtime   !time cloud droplets spend in intermediate 2/3 of cloud
            ELSEIF (l_q(l_i,2) .LE. l_y_base)THEN
               l_q(l_i,6)=C_0_R    !cloud droplets loose memory when "leaving" cloud
               l_hq(l_i,6)=C_0_R   !cloud droplets loose memory when "leaving" cloud
            ENDIF
         ENDDO
      ENDIF


  RETURN
END SUBROUTINE PARTICLE_TIME_RESIDENCE
