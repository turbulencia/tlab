#include "types.h"
#include "avgij_map.h"
      
SUBROUTINE DNS_SAVE_SCBDGIJ_UV1(NNstat, m_u_x, m_v_y, &
     m_rho, m_u, m_v, m_w, m_z1, z1, tmp1, tmp2, tmp3, tmp4, &
     tmp5, tmp6, mean1d_sc, wrk2d)
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL m_u_x(*)
  TREAL m_v_y(*)
  TREAL m_rho(*)
  TREAL m_u(*)
  TREAL m_v(*)
  TREAL m_w(*)
  TREAL m_z1(*)
  TREAL tmp1(*)
  TREAL tmp2(*)
  TREAL tmp3(*)
  TREAL tmp4(*)
  TREAL tmp5(*)
  TREAL tmp6(*)

  TREAL z1(imax,jmax,kmax)
  TREAL mean1d_sc(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)
  TINTEGER j

  CALL REDUCE( imax, jmax, kmax, z1, nstatavg, statavg, m_z1)

#define m_rho_z1_u_x   tmp1
#define m_rho_z1_v_y   tmp2
#define m_rho_z1_2_u_x tmp3
#define m_rho_z1_2_v_y tmp4

  DO j = 1,NNstat*kmax
     m_rho_z1_u_x(j) = m_rho(j)*m_z1(j)*m_u_x(j)
     m_rho_z1_v_y(j) = m_rho(j)*m_z1(j)*m_v_y(j)
     m_rho_z1_2_u_x(j) = m_rho_z1_u_x(j)*m_z1(j)
     m_rho_z1_2_v_y(j) = m_rho_z1_v_y(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_2_u_x, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_2_v_y, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSUx(j) = MS_RSUx(j) + wrk2d(j,1)
     MS_RSVy(j) = MS_RSVy(j) + wrk2d(j,2)
     MS_RSSUx(j) = MS_RSSUx(j) + wrk2d(j,3)
     MS_RSSVy(j) = MS_RSSVy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_z1_u_x
#undef m_rho_z1_v_y
#undef m_rho_z1_2_u_x
#undef m_rho_z1_2_v_y

#define m_rho_u_u_x_z1 tmp1
#define m_rho_v_v_y_z1 tmp2
#define m_rho_u_x_v_z1 tmp3
#define m_rho_u_v_y_z1 tmp4
#define m_rho_u_x_w_z1 tmp5
#define m_rho_v_y_w_z1 tmp6

  DO j = 1,NNstat*kmax
     m_rho_u_u_x_z1(j) = m_rho(j)*m_z1(j)*m_u(j)*m_u_x(j)
     m_rho_v_v_y_z1(j) = m_rho(j)*m_z1(j)*m_v(j)*m_v_y(j)
     m_rho_u_x_v_z1(j) = m_rho(j)*m_z1(j)*m_u_x(j)*m_v(j)
     m_rho_u_v_y_z1(j) = m_rho(j)*m_z1(j)*m_u(j)*m_v_y(j)
     m_rho_u_x_w_z1(j) = m_rho(j)*m_z1(j)*m_u_x(j)*m_w(j)
     m_rho_v_y_w_z1(j) = m_rho(j)*m_z1(j)*m_v_y(j)*m_w(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_u_x_z1, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_v_v_y_z1, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_x_v_z1, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_v_y_z1, wrk2d(1,4), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_u_x_w_z1, wrk2d(1,5), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_v_y_w_z1, wrk2d(1,6), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RUUSx(j) = MS_RUUSx(j) + C_2_R*wrk2d(j,1)
     MS_RVVSy(j) = MS_RVVSy(j) + C_2_R*wrk2d(j,2)
     MS_RUVSx(j) = MS_RUVSx(j) + wrk2d(j,3)
     MS_RUVSy(j) = MS_RUVSy(j) + wrk2d(j,4)
     MS_RUWSx(j) = MS_RUWSx(j) + wrk2d(j,5)
     MS_RVWSy(j) = MS_RVWSy(j) + wrk2d(j,6)
  ENDDO

#undef m_rho_u_u_x_z1
#undef m_rho_v_v_y_z1
#undef m_rho_u_x_v_z1
#undef m_rho_u_v_y_z1
#undef m_rho_u_x_w_z1
#undef m_rho_v_y_w_z1

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_UV1
