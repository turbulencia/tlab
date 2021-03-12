#include "types.h"
#include "avgij_map.h"
      
SUBROUTINE DNS_SAVE_SCBDGIJ_R1(NNstat, m_rho_x, m_rho_y, &
     m_u, m_v, m_w, m_z1, z1, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, mean1d_sc, wrk2d)
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL m_rho_x(*)
  TREAL m_rho_y(*)
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

#define m_rho_x_z1    tmp4
#define m_rho_y_z1    tmp5
#define m_rho_x_z1_z1 tmp2
#define m_rho_y_z1_z1 tmp3

  DO j = 1,NNstat*kmax
     m_rho_x_z1(j) = m_rho_x(j)*m_z1(j)
     m_rho_y_z1(j) = m_rho_y(j)*m_z1(j)
     m_rho_x_z1_z1(j) = m_rho_x_z1(j)*m_z1(j)
     m_rho_y_z1_z1(j) = m_rho_y_z1(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_z1, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_z1, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_SRx(j) = MS_SRx(j) + wrk2d(j,1)
     MS_SRy(j) = MS_SRy(j) + wrk2d(j,2)
     MS_RSSx(j) = MS_RSSx(j) + wrk2d(j,3)
     MS_RSSy(j) = MS_RSSy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_x_z1
#undef m_rho_y_z1
#undef m_rho_x_z1_z1
#undef m_rho_y_z1_z1

#define m_rho_x_z1_u tmp2
#define m_rho_y_z1_v tmp3
#define m_rho_x_z1_2_u tmp4
#define m_rho_y_z1_2_v tmp5

  DO j = 1,NNstat*kmax
     m_rho_x_z1_u(j) = m_rho_x(j)*m_z1(j)*m_u(j)
     m_rho_y_z1_v(j) = m_rho_y(j)*m_z1(j)*m_v(j)
     m_rho_x_z1_2_u(j) = m_rho_x_z1_u(j)*m_z1(j)
     m_rho_y_z1_2_v(j) = m_rho_y_z1_v(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_u, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_v, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_2_u, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_2_v, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSUx(j) = MS_RSUx(j) + wrk2d(j,1)
     MS_RSVy(j) = MS_RSVy(j) + wrk2d(j,2)
     MS_RSSUx(j) = MS_RSSUx(j) + wrk2d(j,3)
     MS_RSSVy(j) = MS_RSSVy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_x_z1_u
#undef m_rho_y_z1_v
#undef m_rho_x_z1_2_u
#undef m_rho_y_z1_2_v

#define m_rho_x_z1_v tmp2
#define m_rho_y_z1_u tmp3
#define m_rho_x_z1_w tmp4
#define m_rho_y_z1_w tmp5

  DO j = 1,NNstat*kmax
     m_rho_x_z1_v(j) = m_rho_x(j)*m_z1(j)*m_v(j)
     m_rho_y_z1_u(j) = m_rho_y(j)*m_z1(j)*m_u(j)
     m_rho_x_z1_w(j) = m_rho_x(j)*m_z1(j)*m_w(j)
     m_rho_y_z1_w(j) = m_rho_y(j)*m_z1(j)*m_w(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_v, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_u, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_w, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_w, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSVx(j) = MS_RSVx(j) + wrk2d(j,1)
     MS_RSUy(j) = MS_RSUy(j) + wrk2d(j,2)
     MS_RSWx(j) = MS_RSWx(j) + wrk2d(j,3)
     MS_RSWy(j) = MS_RSWy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_x_z1_v
#undef m_rho_y_z1_u
#undef m_rho_x_z1_w
#undef m_rho_y_z1_w

#define m_rho_x_u_u_z1 tmp1
#define m_rho_y_v_v_z1 tmp6
#define m_rho_x_u_v_z1 tmp2
#define m_rho_y_u_v_z1 tmp3
#define m_rho_x_u_w_z1 tmp4
#define m_rho_y_v_w_z1 tmp5

  DO j = 1,NNstat*kmax
     m_rho_x_u_u_z1(j) = m_rho_x(j)*m_z1(j)*m_u(j)**2
     m_rho_y_v_v_z1(j) = m_rho_y(j)*m_z1(j)*m_v(j)**2
     m_rho_x_u_v_z1(j) = m_rho_x(j)*m_z1(j)*m_u(j)*m_v(j)
     m_rho_y_u_v_z1(j) = m_rho_y(j)*m_z1(j)*m_u(j)*m_v(j)
     m_rho_x_u_w_z1(j) = m_rho_x(j)*m_z1(j)*m_u(j)*m_w(j)
     m_rho_y_v_w_z1(j) = m_rho_y(j)*m_z1(j)*m_v(j)*m_w(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_u_u_z1, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_v_v_z1, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_u_v_z1, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_u_v_z1, wrk2d(1,4), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_x_u_w_z1, wrk2d(1,5), &
       wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_rho_y_v_w_z1, wrk2d(1,6), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RUUSx(j) = MS_RUUSx(j) + wrk2d(j,1)
     MS_RVVSy(j) = MS_RVVSy(j) + wrk2d(j,2)
     MS_RUVSx(j) = MS_RUVSx(j) + wrk2d(j,3)
     MS_RUVSy(j) = MS_RUVSy(j) + wrk2d(j,4)
     MS_RUWSx(j) = MS_RUWSx(j) + wrk2d(j,5)
     MS_RVWSy(j) = MS_RVWSy(j) + wrk2d(j,6)
  ENDDO

#undef m_rho_x_u_u_z1 
#undef m_rho_y_v_v_z1 
#undef m_rho_x_u_v_z1
#undef m_rho_y_u_v_z1 
#undef m_rho_x_u_w_z1 
#undef m_rho_y_v_w_z1 

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_R1

