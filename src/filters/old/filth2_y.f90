#include "types.h"
#include "dns_error.h"

SUBROUTINE FILTH2_Y(iunify, j1bc, imax, jmax, kmax, ny0, ny1, cf2y, z1, zf1, wrk)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER iunify, j1bc
  TINTEGER imax, jmax, kmax
  TINTEGER ny0, ny1
  TREAL cf2y(*)
  TREAL z1(imax, jmax, kmax)
  TREAL zf1(imax, jmax, kmax)
  TREAL wrk(imax, jmax, kmax)

  TINTEGER nij, nik
  TINTEGER i2

  i2 = 2
  nij = imax*jmax
  nik = imax*kmax

  IF ( MOD(ny0,i2) .NE. 0 .AND. MOD(ny1,i2) .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'FILTH2_Y. NY2 is not even')
     CALL DNS_STOP(DNS_ERROR_LESEVEN)
  ENDIF

  CALL DNS_TRANSPOSE(z1, nij, kmax, nij, zf1, kmax)

  IF ( iunify .EQ. 0 ) THEN
     IF ( j1bc .EQ. 0 ) THEN
        IF ( ny0+ny1 .EQ. 4 ) THEN
           CALL FILTH2MPPD4(jmax, nik, cf2y, zf1, wrk)
        ELSE IF ( ny0+ny1 .EQ. 6 ) THEN
           CALL FILTH2MPPD6(jmax, nik, cf2y, zf1, wrk)
        ELSE
           CALL FILTH2MPPD(jmax, nik, ny0, ny1, cf2y, zf1, wrk)
        ENDIF
     ELSE
        CALL FILTH2MP(jmax, nik, ny0, ny1, cf2y, zf1, wrk)
     ENDIF
  ELSE
     IF ( j1bc .EQ. 0 ) THEN
        CALL FILTH2MPPDNU(jmax, nik, ny0, ny1, cf2y, zf1, wrk)
     ELSE
        CALL FILTH2MPNU(jmax, nik, ny0, ny1, cf2y, zf1, wrk)
     ENDIF
  ENDIF

  CALL DNS_TRANSPOSE(wrk, kmax, nij, kmax, zf1, nij)

  RETURN
END SUBROUTINE FILTH2_Y

