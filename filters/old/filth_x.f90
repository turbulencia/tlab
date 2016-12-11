#include "types.h"
#include "dns_error.h"

SUBROUTINE FILTH_X(iunifx, i1bc, imax, jmax, kmax, nx, cfx, z1, zf1, wrk)
  
  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  TINTEGER iunifx, i1bc
  TINTEGER imax, jmax, kmax
  TINTEGER nx
  TREAL cfx(*)
  TREAL z1(imax, jmax, kmax)
  TREAL zf1(imax, jmax, kmax)
  TREAL wrk(imax, jmax, kmax)

  TINTEGER njk
  TINTEGER i2

  i2 = 2
  njk = jmax*kmax

  IF ( MOD(nx,i2) .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'FILTH_X. NX is not even')
     CALL DNS_STOP(DNS_ERROR_LESEVEN)
  ENDIF

  CALL DNS_TRANSPOSE(z1, imax, njk, imax, zf1, njk)

  IF ( iunifx .EQ. 0 ) THEN
     IF ( i1bc .EQ. 0 ) THEN
        IF ( nx .EQ. 2 ) THEN
           CALL FILTHMPPD2(imax, njk, zf1, wrk)
        ELSE IF ( nx .EQ. 4 ) THEN
           CALL FILTHMPPD4(imax, njk, zf1, wrk)
        ELSE
           CALL FILTHMPPD(imax, njk, nx, zf1, wrk)
        ENDIF
     ELSE
        CALL FILTHMP(imax, njk, nx, cfx, zf1, wrk)
     ENDIF
  ELSE
     IF ( i1bc .EQ. 0 ) THEN
        CALL FILTHMPPDNU(imax, njk, nx, cfx, zf1, wrk)
     ELSE
        IF ( nx .EQ. 2 ) THEN
           CALL FILTHMPNU2(imax, njk, cfx, zf1, wrk)
        ELSE IF ( nx .EQ. 4 ) THEN
           CALL FILTHMPNU4(imax, njk, cfx, zf1, wrk)
        ELSE IF ( nx .EQ. 6 ) THEN
           CALL FILTHMPNU6(imax, njk, cfx, zf1, wrk)
        ELSE
           CALL FILTHMPNU(imax, njk, nx, cfx, zf1, wrk)
        ENDIF
     ENDIF
  ENDIF

  CALL DNS_TRANSPOSE(wrk, njk, imax, njk, zf1, imax)

  RETURN
END SUBROUTINE FILTH_X

