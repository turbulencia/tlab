#include "types.h"
#include "avgij_map.h"

SUBROUTINE DNS_SAVE_SCBDGIJ_UV2(NNstat, m_v_x, m_u_y, &
     m_rho, m_u, m_v, m_w, m_z1, z1, tmp1, tmp2, tmp3, tmp4, &
     mean1d_sc, wrk2d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL m_v_x(*)
  TREAL m_u_y(*)
  TREAL m_rho(*)
  TREAL m_u(*)
  TREAL m_v(*)
  TREAL m_w(*)
  TREAL m_z1(*)
  TREAL tmp1(*)
  TREAL tmp2(*)
  TREAL tmp3(*)
  TREAL tmp4(*)
  TREAL z1(imax,jmax,kmax)
  TREAL mean1d_sc(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)
  TINTEGER j

  CALL REDUCE( imax, jmax, kmax, z1, nstatavg, statavg, m_z1)

#define m_rho_z1_v_x tmp1
#define m_rho_z1_u_y tmp2
#define m_rho_u_v_x_z1 tmp3
#define m_rho_u_y_v_z1 tmp4

  DO j = 1,NNstat*kmax
     m_rho_z1_v_x(j) = m_rho(j)*m_z1(j)*m_v_x(j)
     m_rho_z1_u_y(j) = m_rho(j)*m_z1(j)*m_u_y(j)
     m_rho_u_v_x_z1(j) = m_rho(j)*m_z1(j)*m_u(j)*m_v_x(j)
     m_rho_u_y_v_z1(j) = m_rho(j)*m_z1(j)*m_u_y(j)*m_v(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_v_x_z1, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_y_v_z1, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSVx(j) = MS_RSVx(j) + wrk2d(j,1)
     MS_RSUy(j) = MS_RSUy(j) + wrk2d(j,2)
     MS_RUVSx(j) = MS_RUVSx(j) + wrk2d(j,3)
     MS_RUVSy(j) = MS_RUVSy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_z1_v_x
#undef m_rho_z1_u_y
#undef m_rho_u_v_x_z1
#undef m_rho_u_y_v_z1

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_UV2
