!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/03/25 - J.P. Mellado
!#              Adding possible dependence perpendicular to boundary
!# 2007/10/26 - J.P. Mellado
!#              Jet case
!#
!########################################################################
!# DESCRIPTION
!#
!# The jet case have to be finished. I would do them based on
!# the initial field by taking averages.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE BOUNDARY_INIT_HB(dx,dz, q,s, txc, buffer_hb)

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax
  USE DNS_GLOBAL, ONLY : imode_sim, imode_eqns, icalc_scal
  USE DNS_GLOBAL, ONLY : inb_flow, inb_scal, scalez, area
  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*)                         :: dx, dz
  TREAL, DIMENSION(imax*jmax*kmax,*)          :: q, s, txc
  TREAL, DIMENSION(imax,buff_nps_jmin,kmax,*) :: buffer_hb

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG_IK, AVG1V1D, COV2V2D, COV2V1D
  TINTEGER i, j, is, iq

  TREAL, DIMENSION(:), POINTER :: r_loc, e_loc

! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN; e_loc => txc(:,2); r_loc => q(:,5)
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN; e_loc => q(:,4);   r_loc => q(:,5);
  ELSE;                                                                  r_loc => txc(:,1); ENDIF

! ###################################################################
! Temporally evolving shear layer
! ###################################################################
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN

! Flow variables
     DO iq = 1,3
     DO j = 1, buff_nps_u_jmin
        IF ( buff_hard_on(iq) .EQ. 0 ) &
             buff_hard(iq,2) = COV2V2D(imax,jmax,kmax, j, i1,i1, r_loc,q(1,iq), dx,dz, area)
        buffer_hb(:,j,:,iq) = buff_hard(iq,2)
     ENDDO
     ENDDO

! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1, buff_nps_e_jmin
        IF ( buff_hard_on(4) .EQ. 0 ) &
             buff_hard(4,2) = COV2V2D(imax,jmax,kmax, j, i1,i1, r_loc,e_loc, dx,dz, area)
        buffer_hb(:,j,:,4) = buff_hard(4,2)
        IF ( buff_hard_on(5) .EQ. 0 ) &
             buff_hard(5,2) = AVG_IK(imax,jmax,kmax, j, r_loc, dx,dz, area)
        buffer_hb(:,j,:,5) = buff_hard(5,2)
     ENDDO
     ENDIF

! Scalars
     IF ( icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
     DO j = 1,buff_nps_u_jmin
        IF ( buff_hard_on(inb_flow+is) .EQ. 0 ) &
             buff_hard(inb_flow+is,2) = COV2V2D(imax,jmax,kmax, j, i1,i1, r_loc,s(1,is), dx,dz, area)
        buffer_hb(:,j,:,inb_flow+is) = buff_hard(inb_flow+is,2)
     ENDDO
     ENDDO
     ENDIF

! ###################################################################
! Spatially evolving jet
! ###################################################################
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN  

! Flow variables
     DO iq = 1,3
     DO j = 1, buff_nps_u_jmin; DO i = 1,imax
        IF ( buff_hard_on(iq) .EQ. 0 ) &
             buff_hard(iq,2) = COV2V1D(imax,jmax,kmax, i,j, r_loc,q(1,iq), dz, scalez)
        buffer_hb(i,j,:,iq) = buff_hard(iq,2)
     ENDDO; ENDDO
     ENDDO

! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1, buff_nps_e_jmin; DO i = 1,imax
        IF ( buff_hard_on(4) .EQ. 0 ) &
             buff_hard(4,2) = COV2V1D(imax,jmax,kmax, i,j, r_loc,e_loc, dz, scalez)
        buffer_hb(i,j,:,4) = buff_hard(4,2)
        IF ( buff_hard_on(5) .EQ. 0 ) &
             buff_hard(5,2) = AVG1V1D(imax,jmax,kmax, i,j, r_loc,       dz, scalez)
        buffer_hb(i,j,:,5) = buff_hard(5,2)
     ENDDO; ENDDO
     ENDIF

! Scalars
     IF ( icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
     DO j = 1,buff_nps_u_jmin; DO i = 1,imax
        IF ( buff_hard_on(inb_flow+is) .EQ. 0 ) &
             buff_hard(inb_flow+is,2) = COV2V1D(imax,jmax,kmax, i,j, r_loc,s(1,is), dz, scalez)
        buffer_hb(i,j,:,inb_flow+is) = buff_hard(inb_flow+is,2)
     ENDDO; ENDDO

     ENDDO
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_INIT_HB
