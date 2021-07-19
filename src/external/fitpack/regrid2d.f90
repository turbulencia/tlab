      SUBROUTINE REGRID2D(imax, jmax, x, y, kx, ky, f, iwrk_size, wrk)

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

      TREAL xb, xe, yb, ye
      TREAL s, fp
      TINTEGER iopt, nxest, nyest, lwrk, kwrk
      TINTEGER mx, my, ip(6)
      TINTEGER nx, ny, ier, i
      TINTEGER iwrk_size
      CHARACTER*32 msg

! ###############################################
! # Initialization of SPLINES subroutines

      iopt = 0
      mx = imax
      my = jmax
      xb = x(1)
      xe = x(imax)
      yb = y(1)
      ye = y(jmax)
      s = C_0_R
      nxest = imax + kx + 1
      nyest = jmax + ky + 1
      lwrk = 4+nxest*(jmax+2*kx+5)+nyest*(2*ky+5)+imax*(kx+1)+&
           jmax*(ky+1)+MAX(jmax,nxest)
      kwrk = 3+mx+my+nxest+nyest

! Working array relative possitions
      ip(1) = 13
      ip(2) = ip(1) + nxest
      ip(3) = ip(2) + nyest
      ip(4) = ip(3) + nxest*nyest
      ip(5) = ip(4) + lwrk
      ip(6) = ip(5) + kwrk

      IF ( iwrk_size .LT. ip(6) ) THEN
         CALL TLAB_WRITE_ASCII(efile, 'REGRID2D Error: Insuficient work space')
         CALL TLAB_STOP(DNS_ERROR_WRKOVERFLW)
      ENDIF

      call regrid(iopt, imax, x, jmax, y, f, xb, xe, &
           yb, ye, kx, ky, s, nxest, nyest, nx, wrk(ip(1)), &
           ny, wrk(ip(2)), wrk(ip(3)), fp, wrk(ip(4)), lwrk, &
           wrk(ip(5)), kwrk, ier)
         
      IF ( ier .NE. 0 .AND. ier .NE. -1 ) THEN
         DO i=1, 32
            msg(i:i) = ' '
         ENDDO
         WRITE(msg,*) 'REGRID ERROR = ',ier
         CALL TLAB_WRITE_ASCII(efile, msg)
         CALL TLAB_STOP(DNS_ERROR_REGRID)
      ENDIF

      DO i=1, 6
         wrk(i) = ip(i) + 0.1
      ENDDO
      
      wrk(7) = nx + 0.1
      wrk(8) = ny + 0.1
      wrk(9) = kx + 0.1
      wrk(10) = ky + 0.1
      wrk(11) = lwrk + 0.1
      wrk(12) = kwrk + 0.1

      RETURN
      END


