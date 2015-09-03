#include "types.h"
#include "avgij_map.h"
      
SUBROUTINE DNS_SAVE_SCBDGIJ_P1(NNstat, m_p_x, m_p_y, m_p_z, &
     m_z1, z1, m_p_x_z1, m_p_y_z1, m_p_z_z1, mean1d_sc, wrk2d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL m_p_x(*)
  TREAL m_p_y(*)
  TREAL m_p_z(*)
  TREAL m_z1(*)
  TREAL m_p_x_z1(*)
  TREAL m_p_y_z1(*)
  TREAL m_p_z_z1(*)
  TREAL z1(imax,jmax,kmax)
  TREAL mean1d_sc(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)

  TINTEGER j

  CALL REDUCE( imax, jmax, kmax, z1, nstatavg, statavg, m_z1)

  DO j = 1,NNstat*kmax
     m_p_x_z1(j) = m_p_x(j)*m_z1(j)
     m_p_y_z1(j) = m_p_y(j)*m_z1(j)
     m_p_z_z1(j) = m_p_z(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_p_x_z1, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_p_y_z1, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_p_z_z1, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_SPx(j) = MS_SPx(j) + wrk2d(j,1)
     MS_SPy(j) = MS_SPy(j) + wrk2d(j,2)
     MS_SPz(j) = MS_SPz(j) + wrk2d(j,3)
  ENDDO

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_P1
