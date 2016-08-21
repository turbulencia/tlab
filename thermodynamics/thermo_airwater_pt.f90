#include "types.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/10/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate liquid content from p, T and water content.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_PT(nx, ny, nz, s, p, T)

  USE THERMO_GLOBAL, ONLY : GRATIO, MRATIO, WGHT_INV, THERMO_AI, dsmooth

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: p,T
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(OUT) :: s

! -------------------------------------------------------------------
  TINTEGER ij
  TREAL qsat, dqldqt, dsmooth_loc

! ###################################################################
  CALL THERMO_POLYNOMIAL_PSAT(nx, ny, nz, T, s(1,2))
  DO ij = 1,nx*ny*nz
! this is really the vapor content
     qsat = C_1_R/(MRATIO*p(ij)/s(ij,2)-C_1_R)*WGHT_INV(2)/WGHT_INV(1)*(1-s(ij,1))
     IF ( qsat .GE. s(ij,1) ) THEN
        s(ij,2) = C_0_R
     ELSE
        s(ij,2) = s(ij,1) - qsat
     ENDIF
     IF ( dsmooth .GT. C_0_R ) THEN
        qsat = qsat/(1-s(ij,1))
        dqldqt = C_1_R + qsat
! this is the real qsat
        qsat = qsat/(C_1_R+qsat)
        dsmooth_loc = dsmooth*qsat
        s(ij,2) = dsmooth_loc*dqldqt&
             *LOG(EXP((s(ij,1)-qsat)/dsmooth_loc)+C_1_R)
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE THERMO_AIRWATER_PT
