SUBROUTINE FILTH_Y(iunify, j1bc, imax, jmax, kmax, ny, cfy, z1, zf1, wrk)
  
#include "types.h"
#include "dns_error.h"

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER iunify, j1bc
  TINTEGER imax, jmax, kmax
  TINTEGER ny
  TREAL cfy(*)
  TREAL z1(imax, jmax, kmax)
  TREAL zf1(imax, jmax, kmax)
  TREAL wrk(imax, jmax, kmax)

  TINTEGER nij, nik
  TINTEGER i2

  i2 = 2
  nij = imax*jmax
  nik = imax*kmax

  IF ( MOD(ny,i2) .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'FILTH_Y. NY is not even')
     CALL DNS_STOP(DNS_ERROR_LESEVEN)
  ENDIF

  CALL DNS_TRANSPOSE(z1, nij, kmax, nij, zf1, kmax)

  IF ( iunify .EQ. 0 ) THEN
     IF ( j1bc .EQ. 0 ) THEN
        IF ( ny .EQ. 2 ) THEN
           CALL FILTHMPPD2(jmax, nik, zf1, wrk)
        ELSE IF ( ny .EQ. 4 ) THEN
           CALL FILTHMPPD4(jmax, nik, zf1, wrk)
        ELSE
           CALL FILTHMPPD(jmax, nik, ny, zf1, wrk)
        ENDIF
     ELSE
        CALL FILTHMP(jmax, nik, ny, cfy, zf1, wrk)
     ENDIF
  ELSE
     IF ( j1bc .EQ. 0 ) THEN
        CALL FILTHMPPDNU(jmax, nik, ny, cfy, zf1, wrk)
     ELSE
        IF ( ny .EQ. 2 ) THEN
           CALL FILTHMPNU2(jmax, nik, cfy, zf1, wrk)
        ELSE IF ( ny .EQ. 4 ) THEN
           CALL FILTHMPNU4(jmax, nik, cfy, zf1, wrk)
        ELSE IF ( ny .EQ. 6 ) THEN
           CALL FILTHMPNU6(jmax, nik, cfy, zf1, wrk)
        ELSE
           CALL FILTHMPNU(jmax, nik, ny, cfy, zf1, wrk)
        ENDIF
     ENDIF
  ENDIF

  CALL DNS_TRANSPOSE(wrk, kmax, nij, kmax, zf1, nij)

  RETURN
END SUBROUTINE FILTH_Y

