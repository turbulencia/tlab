#include "types.h"
#include "dns_error.h"

SUBROUTINE THERMO_POLYNOMIAL_PSAT(nx,ny,nz, T, p)

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
        p(i) = THERMO_PSAT(NPSAT)
        DO ipsat = NPSAT-1,1,-1
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

! ###################################################################
! ###################################################################
SUBROUTINE THERMO_POLYNOMIAL_DPSAT(nx,ny,nz, T, dp)

  USE THERMO_GLOBAL, ONLY : THERMO_PSAT, NPSAT

  IMPLICIT NONE

  TINTEGER nx, ny, nz
  TREAL T(*)
  TREAL dp(*)

! -------------------------------------------------------------------
  TINTEGER i, ipsat

! ###################################################################
  IF ( NPSAT .GT. 0 ) THEN
     DO i = 1, nx*ny*nz
        dp(i) = C_0_R
        DO ipsat = NPSAT-1,1,-1
           dp(i) = dp(i)*T(i) + THERMO_PSAT(ipsat+1) *M_REAL(ipsat)
        ENDDO
     ENDDO
  ELSE
     DO i = 1, nx*ny*nz
        dp(i) = C_0_R
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE THERMO_POLYNOMIAL_DPSAT
