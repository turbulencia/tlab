SUBROUTINE DNS_SAVE_EXTRA_I( icount, ie, evars, p_dat)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "types.h"

  TINTEGER icount
  TREAL evars(imax, jmax, kmax)

  TREAL p_dat(nspa_rest,nstatpln,jmax,nstatplnvars,kmax)

  TINTEGER i, j, k, n, ip, ie

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING DNS_SAVE_EXTRA_I' )
#endif

  ip = inb_vars + 1 + ie 
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = evars(i,j,k)
        ENDDO
     ENDDO
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING DNS_SAVE_EXTRA_I' )
#endif

  RETURN
END SUBROUTINE DNS_SAVE_EXTRA_I
