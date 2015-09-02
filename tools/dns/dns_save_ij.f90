SUBROUTINE DNS_SAVE_IJ(icount, rho, u, v, w, p, z1, l_dat)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "types.h"

  TINTEGER icount
  TREAL, DIMENSION(imax,jmax,kmax)   :: rho, u, v, w, p
  TREAL, DIMENSION(imax,jmax,kmax,*) :: z1

  TREAL l_dat(nspa_rest,nstatlin,inb_vars,kmax)

  TINTEGER i, j, k, n, is, ip

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING DNS_SAVE_IJ' )
#endif

  ip = 0

  ip = ip + 1
  DO k = 1,kmax
     DO n = 1,nstatlin
        i = statlin_i(n)
        j = statlin_j(n)
        l_dat(icount,n,ip,k) = u(i,j,k)
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO n = 1,nstatlin
        i = statlin_i(n)
        j = statlin_j(n)
        l_dat(icount,n,ip,k) = v(i,j,k)
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO n = 1,nstatlin
        i = statlin_i(n)
        j = statlin_j(n)
        l_dat(icount,n,ip,k) = w(i,j,k)
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO n = 1,nstatlin
        i = statlin_i(n)
        j = statlin_j(n)
        l_dat(icount,n,ip,k) = p(i,j,k)
     ENDDO
  ENDDO

  ip = ip + 1
  DO k = 1,kmax
     DO n = 1,nstatlin
        i = statlin_i(n)
        j = statlin_j(n)
        l_dat(icount,n,ip,k) = rho(i,j,k)
     ENDDO
  ENDDO

  DO is=1, inb_scal
     ip = ip + 1
     DO k = 1,kmax
        DO n = 1,nstatlin
           i = statlin_i(n)
           j = statlin_j(n)
           l_dat(icount,n,ip,k) = z1(i,j,k,is)
        ENDDO
     ENDDO
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING DNS_SAVE_IJ' )
#endif

  RETURN
END SUBROUTINE DNS_SAVE_IJ
