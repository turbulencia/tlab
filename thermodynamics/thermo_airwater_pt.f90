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
SUBROUTINE THERMO_AIRWATER_PT(nx, ny, nz, z1, p, T)

  USE THERMO_GLOBAL, ONLY : GRATIO, MRATIO, WGHT_INV, THERMO_AI, dsmooth

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL z1(nx*ny*nz,*), T(*), p(*)

! -------------------------------------------------------------------
  TINTEGER ij
  TREAL qsat, dqldqt, dsmooth_loc

! ###################################################################
  CALL THERMO_POLYNOMIAL_PSAT(nx, ny, nz, T, z1(1,2))
  DO ij = 1,nx*ny*nz
! this is really the vapor content
     qsat = C_1_R/(MRATIO*p(ij)/z1(ij,2)-C_1_R)*WGHT_INV(2)/WGHT_INV(1)*(1-z1(ij,1))
     IF ( qsat .GE. z1(ij,1) ) THEN
        z1(ij,2) = C_0_R
     ELSE
        z1(ij,2) = z1(ij,1) - qsat
     ENDIF
     IF ( dsmooth .GT. C_0_R ) THEN
        qsat = qsat/(1-z1(ij,1))
        dqldqt = C_1_R + qsat
! this is the real qsat
        qsat = qsat/(C_1_R+qsat)
        dsmooth_loc = dsmooth*qsat
        z1(ij,2) = dsmooth_loc*dqldqt&
             *LOG(EXP((z1(ij,1)-qsat)/dsmooth_loc)+C_1_R)
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE THERMO_AIRWATER_PT
