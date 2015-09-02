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

SUBROUTINE BOUNDARY_INIT_HT(dx,dz, q,s, txc, buffer_ht)

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax
  USE DNS_GLOBAL, ONLY : imode_sim, imode_eqns, icalc_scal
  USE DNS_GLOBAL, ONLY : inb_flow, inb_scal, scalez, area
  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*)                         :: dx, dz
  TREAL, DIMENSION(imax*jmax*kmax,*)          :: q, s, txc
  TREAL, DIMENSION(imax,buff_nps_jmax,kmax,*) :: buffer_ht

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG_IK, AVG1V1D, AVG2V2D, COV2V1D
  TINTEGER i, j, is, iq, jloc

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
     DO j = 1, buff_nps_u_jmax; jloc = buff_u_jmax + j - 1
        IF ( buff_hard_on(iq) .EQ. 0 ) &
             buff_hard(iq,1) = AVG2V2D(imax,jmax,kmax, jloc, r_loc,q(1,iq))
        buffer_ht(:,j,:,iq) = buff_hard(iq,1)
     ENDDO
     ENDDO

! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1, buff_nps_e_jmax; jloc = buff_e_jmax + j - 1
        IF ( buff_hard_on(4) .EQ. 0 ) &
             buff_hard(4,1) = AVG2V2D(imax,jmax,kmax, jloc, r_loc,e_loc)
        buffer_ht(:,j,:,4) = buff_hard(4,1)
        IF ( buff_hard_on(5) .EQ. 0 ) &
             buff_hard(5,1) = AVG_IK(imax,jmax,kmax, jloc, r_loc, dx,dz, area)
        buffer_ht(:,j,:,5) = buff_hard(5,1)
     ENDDO
     ENDIF

! Scalars
     IF ( icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
     DO j = 1,buff_nps_u_jmax; jloc = buff_u_jmax + j - 1
        IF ( buff_hard_on(inb_flow+is) .EQ. 0 ) &
             buff_hard(inb_flow+is,1) = AVG2V2D(imax,jmax,kmax, jloc, r_loc,s(1,is))
           buffer_ht(:,j,:,inb_flow+is) = buff_hard(inb_flow+is,1)
     ENDDO
     ENDDO
     ENDIF

! ###################################################################
! Spatially evolving jet
! ###################################################################
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN  

! Flow variables
     DO iq = 1,3
     DO j = 1, buff_nps_u_jmax; jloc = buff_u_jmax + j - 1; DO i = 1,imax
        IF ( buff_hard_on(iq) .EQ. 0 ) &
             buff_hard(iq,1) = COV2V1D(imax,jmax,kmax, i,jloc, r_loc,q(1,iq))
        buffer_ht(i,j,:,iq) = buff_hard(iq,1)
     ENDDO; ENDDO
     ENDDO

! if compressible; energy can have a different buffer extent
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1, buff_nps_e_jmax; jloc = buff_e_jmax + j - 1; DO i = 1,imax
        IF ( buff_hard_on(4) .EQ. 0 ) &
             buff_hard(4,1) = COV2V1D(imax,jmax,kmax, i,jloc, r_loc,e_loc)
        buffer_ht(i,j,:,4) = buff_hard(4,1)
        IF ( buff_hard_on(5) .EQ. 0 ) &
             buff_hard(5,1) = AVG1V1D(imax,jmax,kmax, i,jloc, i1, r_loc)
        buffer_ht(i,j,:,5) = buff_hard(5,1)
     ENDDO; ENDDO
     ENDIF

! Scalars
     IF ( icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
     DO j = 1,buff_nps_u_jmax; jloc = buff_u_jmax + j - 1; DO i = 1,imax
        IF ( buff_hard_on(inb_flow+is) .EQ. 0 ) &
             buff_hard(inb_flow+is,1) = COV2V1D(imax,jmax,kmax, i,jloc, r_loc,s(1,is))
        buffer_ht(i,j,:,inb_flow+is) = buff_hard(inb_flow+is,1)
     ENDDO; ENDDO

     ENDDO
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_INIT_HT
