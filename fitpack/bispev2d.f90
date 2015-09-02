      SUBROUTINE BISPEV2D(imax, jmax, x, y, f, wrk)

#include "types.h"
#include "dns_error.h"

        USE DNS_CONSTANTS, ONLY : efile

      IMPLICIT NONE

      TINTEGER imax, jmax
      TINTEGER kx, ky
      TREAL x(imax)
      TREAL y(jmax)
      TREAL f(imax,jmax)
      TREAL wrk(*)

      TINTEGER lwrk, kwrk
      TINTEGER ip(6)
      TINTEGER nx, ny, ier, i
      CHARACTER*32 msg

! Working array relative possitions
      DO i=1, 6
         ip(i) = INT(wrk(i)+0.1)
      ENDDO
      
      nx = INT(wrk(7) + 0.1)
      ny = INT(wrk(8) + 0.1)
      kx = INT(wrk(9) + 0.1)
      ky = INT(wrk(10) + 0.1)
      lwrk = INT(wrk(11) + 0.1)
      kwrk = INT(wrk(12) + 0.1)
     
      call bispev(wrk(ip(1)), nx, wrk(ip(2)), ny, wrk(ip(3)), kx, ky,&
           x,imax,y,jmax,f,wrk(ip(4)),lwrk,wrk(ip(5)), kwrk, ier)

      IF ( ier .NE. 0 .AND. ier .NE. -1 ) THEN
         DO i=1, 32
            msg(i:i) = ' '
         ENDDO
         WRITE(msg,*) 'BISPEV ERROR = ',ier
         CALL IO_WRITE_ASCII(efile, msg)
         CALL DNS_STOP(DNS_ERROR_REGRID)
      ENDIF

      RETURN
      END


