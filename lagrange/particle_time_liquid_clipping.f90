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
SUBROUTINE PARTICLE_TIME_LIQUID_CLIPPING(s, wrk2d,wrk3d, l_txc, l_tags, l_hq, l_q, l_comm)

  USE DNS_TYPES,  ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, inb_particle, isize_particle
  USE DNS_GLOBAL, ONLY : isize_field, inb_scal_array
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_size_p, ims_pro
#endif

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), TARGET :: s
  TREAL, DIMENSION(*)                     :: wrk2d, wrk3d

  TREAL, DIMENSION(isize_particle,inb_particle) :: l_q, l_hq
  TREAL, DIMENSION(isize_particle), TARGET      :: l_txc
  TREAL, DIMENSION(*)                           :: l_comm
  INTEGER(8), DIMENSION(*)                      :: l_tags
  TINTEGER is, l_i, particle_number_local

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(1) :: data
  TYPE(pointers_dt),   DIMENSION(1) :: data_out

#ifdef USE_MPI
  particle_number_local = ims_size_p(ims_pro+1)
#else
  particle_number_local = INT(particle_number)
#endif

! ###################################################################
! IF negative liquid set lagrange liquid 0
! ###################################################################
  DO is=4,inb_particle_evolution
     DO l_i=1,particle_number_local
        IF( l_q(l_i,is) .LT. 0 ) THEN
           l_q(l_i,is)=C_0_R
           l_q(l_i,6)=C_0_R
           l_hq(l_i,6)=C_0_R
        ENDIF
     ENDDO
  ENDDO

! ###################################################################
! If no liquid around in Eulerian, set liquid droplet to zero
! ###################################################################
  nvar = 0
  nvar = nvar+1; data(nvar)%field(1:imax,1:jmax,1:kmax) => s(:,inb_scal_array); data_out(nvar)%field => l_txc(:)
  CALL FIELD_TO_PARTICLE(nvar, data, data_out, l_q,l_tags,l_comm, wrk2d,wrk3d)

  DO l_i=1,particle_number_local
     IF (l_txc(l_i) .LT. 0.00001) THEN
        DO is=4,inb_particle_evolution
           l_q(l_i,is)=C_0_R
           l_q(l_i,6)=C_0_R
           l_hq(l_i,6)=C_0_R
        ENDDO
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE PARTICLE_TIME_LIQUID_CLIPPING
