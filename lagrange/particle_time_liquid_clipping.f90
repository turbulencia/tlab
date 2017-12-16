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
!# Sets negative particle liquid to zero
!# Sets particle liquid with no eulerian liquid surrounded to zero 
!#
!########################################################################
!# ARGUMENTS 
!#
!# 
!#
!########################################################################
SUBROUTINE PARTICLE_TIME_LIQUID_CLIPPING(s,wrk1d,wrk2d,wrk3d, l_txc, l_tags, l_hq, l_q)    

  USE DNS_GLOBAL, ONLY : inb_particle, isize_particle
  USE DNS_GLOBAL, ONLY : isize_field, inb_scal_array
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_size_p, ims_pro
#endif

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*) :: s
  TREAL, DIMENSION(*)             :: wrk1d, wrk2d, wrk3d

  TREAL, DIMENSION(isize_particle,inb_particle) :: l_q, l_hq
  TREAL, DIMENSION(*)                :: l_txc
  INTEGER(8), DIMENSION(*)           :: l_tags
  TINTEGER is, l_i, local_isize_particle

#ifdef USE_MPI
   local_isize_particle = ims_size_p(ims_pro+1)
#else
   local_isize_particle = particle_number
#endif


! ###################################################################
! IF negative liquid set lagrange liquid 0
! ###################################################################
      IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
         DO is=4,inb_particle_evolution
            DO l_i=1,local_isize_particle
               IF( l_q(l_i,is) .LT. 0 ) THEN
                  l_q(l_i,is)=C_0_R
                  l_q(l_i,6)=C_0_R
                  l_hq(l_i,6)=C_0_R
               ENDIF
            ENDDO
         ENDDO

      ELSE
         DO is=4,inb_particle_evolution
            DO l_i=1,local_isize_particle
               IF( l_q(l_i,is) .LT. 0 ) THEN
                  l_q(l_i,is)=C_0_R
               ENDIF
            ENDDO
         ENDDO

      ENDIF

! ###################################################################
! If no liquid around in Eulerian, set liquid droplet to zero
! ###################################################################
      CALL FIELD_TO_PARTICLE (s(1,inb_scal_array),wrk1d,wrk2d,wrk3d, l_txc, l_tags, l_hq, l_q)  !Update the liquid function
      IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
         DO l_i=1,local_isize_particle
            IF (l_txc(l_i) .LT. 0.00001) THEN
               DO is=4,inb_particle_evolution
                  l_q(l_i,is)=C_0_R
                  l_q(l_i,6)=C_0_R
                  l_hq(l_i,6)=C_0_R
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO l_i=1,local_isize_particle
            IF (l_txc(l_i) .LT. 0.00001) THEN
               DO is=4,inb_particle_evolution
                  l_q(l_i,is)=C_0_R
               ENDDO
            ENDIF
         ENDDO
      ENDIF
 

  RETURN
END SUBROUTINE PARTICLE_TIME_LIQUID_CLIPPING
