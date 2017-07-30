#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Setting buffer zone in output boundary. It only makes sense in
!# spatially evolving cases.
!#
!########################################################################
SUBROUTINE BOUNDARY_INIT_VO(q,s, txc, buffer_q, buffer_s)

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax
  USE DNS_GLOBAL, ONLY : imode_eqns, inb_scal, inb_vars
  USE DNS_LOCAL,  ONLY : BuffFlowImax, BuffScalImax

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax*jmax*kmax,*), INTENT(IN)  :: q, s, txc
  TREAL, DIMENSION(BuffFlowImax%size,jmax,kmax,*) :: buffer_q
  TREAL, DIMENSION(BuffScalImax%size,jmax,kmax,*) :: buffer_s

  TARGET txc, q

! -------------------------------------------------------------------
  TREAL AVG1V1D, COV2V1D
  TREAL dbuff(inb_vars)
  TINTEGER i, j, iq, is, iloc

  TREAL, DIMENSION(:), POINTER :: r_loc, e_loc

! ###################################################################
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN; e_loc => txc(:,2); r_loc => q(:,5)
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN; e_loc => q(:,4);   r_loc => q(:,5);
  ELSE;                                                                  r_loc => txc(:,1); ENDIF

! ###################################################################
! Shear layer and jet profile
! ###################################################################
  DO iq = 1,3
     DO j = 1,jmax
        DO i = 1,BuffFlowImax%size
           iloc = imax -BuffScalImax%size +i
           dbuff(iq) = COV2V1D(imax,jmax,kmax, iloc,j, r_loc,q(1,iq))
           buffer_q(i,j,:,iq) = dbuff(iq)
        ENDDO
     ENDDO
  ENDDO
  
! if compressible
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     DO j = 1,jmax
        DO i = 1,BuffFlowImax%size
           iloc = imax -BuffScalImax%size +i
           dbuff(4) = COV2V1D(imax,jmax,kmax, iloc,j, r_loc,e_loc)
           buffer_q(i,j,:,4) = dbuff(4)
           dbuff(5) = AVG1V1D(imax,jmax,kmax, iloc,j, i1, r_loc)
           buffer_q(i,j,:,5) = dbuff(5)
        ENDDO
     ENDDO
  ENDIF
  
  IF ( BuffScalImax%size .GT. 0 ) THEN
     DO is = 1,inb_scal
        DO j = 1,jmax
           DO i = 1,BuffScalImax%size
              iloc = imax -BuffScalImax%size +i
              dbuff(is) = COV2V1D(imax,jmax,kmax, iloc,j, r_loc,s(1,is))
              buffer_s(i,j,:,is) = dbuff(is)
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE BOUNDARY_INIT_VO

