#include "types.h"
#include "dns_error.h"

SUBROUTINE FILTH2_X(iunifx, i1bc, imax, jmax, kmax, nx0, nx1, cf2x, z1, zf1, wrk)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER iunifx, i1bc
  TINTEGER imax, jmax, kmax
  TINTEGER nx0, nx1
  TREAL cf2x(*)
  TREAL z1(imax, jmax, kmax)
  TREAL zf1(imax, jmax, kmax)
  TREAL wrk(imax, jmax, kmax)

  TINTEGER njk
  TINTEGER i2

  i2 = 2
  njk = jmax*kmax

  IF ( MOD(nx0,i2) .NE. 0 .AND. MOD(nx1,i2) .NE. 0 ) THEN
     CALL TLAB_WRITE_ASCII(efile, 'FILTH2_X. NX2 is not even')
     CALL DNS_STOP(DNS_ERROR_LESEVEN)
  ENDIF

  CALL DNS_TRANSPOSE(z1, imax, njk, imax, zf1, njk)

  IF ( iunifx .EQ. 0 ) THEN
     IF ( i1bc .EQ. 0 ) THEN
        IF ( nx0+nx1 .EQ. 4 ) THEN
           CALL FILTH2MPPD4(imax, njk, cf2x, zf1, wrk)
        ELSE IF ( nx0+nx1 .EQ. 6 ) THEN
           CALL FILTH2MPPD6(imax, njk, cf2x, zf1, wrk)
        ELSE
           CALL FILTH2MPPD(imax, njk, nx0, nx1, cf2x, zf1, wrk)
        ENDIF
     ELSE
        CALL FILTH2MP(imax, njk, nx0, nx1, cf2x, zf1, wrk)
     ENDIF
  ELSE
     IF ( i1bc .EQ. 0 ) THEN
        CALL FILTH2MPPDNU(imax, njk, nx0, nx1, cf2x, zf1, wrk)
     ELSE
        CALL FILTH2MPNU(imax, njk, nx0, nx1, cf2x, zf1, wrk)
     ENDIF
  ENDIF

  CALL DNS_TRANSPOSE(wrk, njk, imax, njk, zf1, imax)

  RETURN
END SUBROUTINE FILTH2_X

