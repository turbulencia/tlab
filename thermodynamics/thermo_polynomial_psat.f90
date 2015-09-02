#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_POLYNOMIAL_PSAT(nx, ny, nz, T, p)

  USE THERMO_GLOBAL, ONLY : THERMO_PSAT, NPSAT

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL T(*)
  TREAL p(*)

! -------------------------------------------------------------------
  TINTEGER i, ipsat

! ###################################################################
  IF ( NPSAT .GT. 0 ) THEN
     DO i = 1, nx*ny*nz
        p(i) = C_0_R
        DO ipsat = NPSAT,1,-1
           p(i) = p(i)*T(i) + THERMO_PSAT(ipsat)
        ENDDO
     ENDDO
  ELSE
     DO i = 1, nx*ny*nz
        p(i) = C_0_R
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE THERMO_POLYNOMIAL_PSAT
