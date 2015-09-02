SUBROUTINE DNS_SAVE_I(icount, rho, u, v, w, p, z1, T, p_dat)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "types.h"

  TINTEGER icount
  TREAL, DIMENSION(imax,jmax,kmax)   :: rho, u, v, w, p, T
  TREAL, DIMENSION(imax,jmax,kmax,*) :: z1

  TREAL p_dat(nspa_rest,nstatpln,jmax,nstatplnvars,kmax)

  TINTEGER i, j, k, n, is, ip

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING DNS_SAVE_I' )
#endif

  ip = 0

  ip = ip + 1
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = u(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = v(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = w(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = p(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = rho(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  DO is=1, inb_scal
     ip = ip + 1
     DO k = 1,kmax
        DO j = 1,jmax
           DO n = 1,nstatpln
              i = statpln(n)
              p_dat(icount,n,j,ip,k) = z1(i,j,k,is)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO j = 1,jmax
        DO n = 1,nstatpln
           i = statpln(n)
           p_dat(icount,n,j,ip,k) = T(i,j,k)
        ENDDO
     ENDDO
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING DNS_SAVE_I' )
#endif

  RETURN
END SUBROUTINE DNS_SAVE_I
