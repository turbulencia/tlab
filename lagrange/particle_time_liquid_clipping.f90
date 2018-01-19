#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

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
SUBROUTINE PARTICLE_TIME_LIQUID_CLIPPING(s, l_q,l_hq,l_txc,l_comm, wrk2d,wrk3d)

  USE DNS_TYPES,  ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_particle
  USE DNS_GLOBAL, ONLY : isize_field, inb_scal_array
  USE LAGRANGE_GLOBAL, ONLY : inb_particle_evolution, particle_number_local, l_g

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), TARGET  :: s

  TREAL, DIMENSION(isize_particle,*)       :: l_q, l_hq
  TREAL, DIMENSION(isize_particle), TARGET :: l_txc
  TREAL, DIMENSION(*)                      :: l_comm
  
  TREAL, DIMENSION(*)                      :: wrk2d, wrk3d

  TINTEGER is, i

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(1) :: data
  TYPE(pointers_dt),   DIMENSION(1) :: data_out

! ###################################################################
! IF negative liquid set lagrange liquid 0
! ###################################################################
  DO is=4,inb_particle_evolution
     DO i=1,particle_number_local
        IF ( l_q(i,is) .LT. C_0_R ) THEN
           l_q(i,is)=C_0_R
           l_q(i,6)=C_0_R
           l_hq(i,6)=C_0_R
        ENDIF
     ENDDO
  ENDDO

! ###################################################################
! If no liquid around in Eulerian, set liquid droplet to zero
! ###################################################################
  nvar = 0
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => s(:,inb_scal_array); data_out(nvar)%field => l_txc(:)
  l_txc = C_0_R
  CALL FIELD_TO_PARTICLE(nvar, data, data_out, l_g,l_q,l_comm, wrk2d,wrk3d)

  DO i=1,particle_number_local
     IF (l_txc(i) .LT. 0.00001) THEN
        DO is=4,inb_particle_evolution
           l_q(i,is)=C_0_R
           l_q(i,6)=C_0_R
           l_hq(i,6)=C_0_R
        ENDDO
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE PARTICLE_TIME_LIQUID_CLIPPING
