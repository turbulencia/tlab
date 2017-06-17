#include "types.h"
#include "avgij_map.h"

SUBROUTINE DNS_SAVE_SCBDGIJ_W1(NNstat, m_w_x, m_w_y, &
     m_rho, m_u, m_v, m_z1, z1, tmp1, tmp2, tmp3, &
     tmp4, mean1d_sc, wrk2d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL m_w_x(*)
  TREAL m_w_y(*)
  TREAL m_rho(*)
  TREAL m_u(*)
  TREAL m_v(*)
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

#define m_rho_z1_w_x   tmp1
#define m_rho_z1_w_y   tmp2
#define m_rho_u_w_x_z1 tmp3
#define m_rho_v_w_y_z1 tmp4

  DO j = 1,NNstat*kmax
     m_rho_z1_w_x(j) = m_rho(j)*m_z1(j)*m_w_x(j)
     m_rho_z1_w_y(j) = m_rho(j)*m_z1(j)*m_w_y(j)
     m_rho_u_w_x_z1(j) = m_rho(j)*m_z1(j)*m_u(j)*m_w_x(j)
     m_rho_v_w_y_z1(j) = m_rho(j)*m_z1(j)*m_v(j)*m_w_y(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_w_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_w_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_w_x_z1, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_v_w_y_z1, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSWx(j) = MS_RSWx(j) + wrk2d(j,1)
     MS_RSWy(j) = MS_RSWy(j) + wrk2d(j,2)
     MS_RUWSx(j) = MS_RUWSx(j) + wrk2d(j,3)
     MS_RVWSy(j) = MS_RVWSy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_z1_w_x
#undef m_rho_z1_w_y
#undef m_rho_u_w_x_z1
#undef m_rho_v_w_y_z1

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_W1
