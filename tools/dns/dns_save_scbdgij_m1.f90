      SUBROUTINE DNS_SAVE_SCBDGIJ_M1(NNstat, m_rho, m_u, m_v, m_w, &
           m_z1, z1, tmp1, tmp2, tmp3, tmp4, tmp5, mean1d_sc, wrk2d)

#include "types.h"
#include "avgij_map.h"
      
  USE DNS_GLOBAL

      IMPLICIT NONE

      TINTEGER NNstat
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
      TREAL z1(imax,jmax,kmax)
      TREAL mean1d_sc(nstatavg,jmax,*)
      TREAL wrk2d(isize_wrk2d,*)
      TINTEGER j

      TREAL inter1, inter2

      CALL REDUCE( imax, jmax, kmax, z1, nstatavg, statavg, m_z1)

#define m_rho_z1    tmp1
#define m_rho_z1_u  tmp2
#define m_rho_z1_v  tmp3
#define m_rho_z1_w  tmp4
#define m_rho_z1_z1 tmp5

      DO j = 1,NNstat*kmax
         m_rho_z1(j) = m_rho(j)*m_z1(j)
         m_rho_z1_u(j) = m_rho_z1(j)*m_u(j)
         m_rho_z1_v(j) = m_rho_z1(j)*m_v(j)
         m_rho_z1_w(j) = m_rho_z1(j)*m_w(j)
         m_rho_z1_z1(j) = m_rho_z1(j)*m_z1(j)
      ENDDO

      CALL SUM_K_V( NNstat, kmax, m_z1, wrk2d(1,1), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1, wrk2d(1,2), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_u, wrk2d(1,3), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_v, wrk2d(1,4), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_w, wrk2d(1,5), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_z1, wrk2d(1,6), &
           wrk2d(1,11))

      DO j = 1,NNstat
         MS_S(j) = MS_S(j) + wrk2d(j,1)
         MS_RS(j) = MS_RS(j) + wrk2d(j,2)
         MS_RSU(j) = MS_RSU(j) + wrk2d(j,3)
         MS_RSV(j) = MS_RSV(j) + wrk2d(j,4)
         MS_RSW(j) = MS_RSW(j) + wrk2d(j,5)
         MS_RSS(j) = MS_RSS(j) + wrk2d(j,6)
      ENDDO

! # Turbulent scalar transport

#define m_rho_z1_u_2 m_rho_z1_u
#define m_rho_z1_v_2 m_rho_z1_v
#define m_rho_z1_w_2 m_rho_z1_w

      DO j = 1,NNstat*kmax
         m_rho_z1_u_2(j) = m_rho_z1_u(j)*m_u(j)
         m_rho_z1_v_2(j) = m_rho_z1_v(j)*m_v(j)
         m_rho_z1_w_2(j) = m_rho_z1_w(j)*m_w(j)
      ENDDO

      CALL SUM_K_V( NNstat, kmax, m_rho_z1_u_2, wrk2d(1,1), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_v_2, wrk2d(1,2), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_w_2, wrk2d(1,3), &
           wrk2d(1,11))

      DO j = 1,NNstat
         MS_RUUS(j) = MS_RUUS(j) + wrk2d(j,1)
         MS_RVVS(j) = MS_RVVS(j) + wrk2d(j,2)
         MS_RWWS(j) = MS_RWWS(j) + wrk2d(j,3)
      ENDDO

#undef m_rho_z1_u_2
#undef m_rho_z1_v_2
#undef m_rho_z1_w_2

#define m_rho_z1_u_v m_rho_z1_u
#define m_rho_z1_u_w m_rho_z1_v
#define m_rho_z1_v_w m_rho_z1_w

      DO j = 1,NNstat*kmax
         m_rho_z1_u_v(j) = m_rho_z1(j)*m_u(j)*m_v(j)
         m_rho_z1_u_w(j) = m_rho_z1(j)*m_u(j)*m_w(j)
         m_rho_z1_v_w(j) = m_rho_z1(j)*m_v(j)*m_w(j)
      ENDDO

      CALL SUM_K_V( NNstat, kmax, m_rho_z1_u_v, wrk2d(1,1), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_u_w, wrk2d(1,2), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_z1_v_w, wrk2d(1,3), &
           wrk2d(1,11))

      DO j = 1,NNstat
         MS_RUVS(j) = MS_RUVS(j) + wrk2d(j,1)
         MS_RUWS(j) = MS_RUWS(j) + wrk2d(j,2)
         MS_RVWS(j) = MS_RVWS(j) + wrk2d(j,3)
      ENDDO

#undef m_rho_z1_u_v
#undef m_rho_z1_u_w
#undef m_rho_z1_v_w

#undef m_rho_z1_u
#undef m_rho_z1_v
#undef m_rho_z1_w

#define m_rho_u_z1_2 tmp2
#define m_rho_v_z1_2 tmp3
#define m_rho_w_z1_2 tmp4

      DO j = 1,NNstat*kmax
         m_rho_u_z1_2(j) = m_rho_z1(j)*m_u(j)*m_z1(j)
         m_rho_v_z1_2(j) = m_rho_z1(j)*m_v(j)*m_z1(j)
         m_rho_w_z1_2(j) = m_rho_z1(j)*m_w(j)*m_z1(j)
      ENDDO

      CALL SUM_K_V( NNstat, kmax, m_rho_u_z1_2, wrk2d(1,1), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_v_z1_2, wrk2d(1,2), &
           wrk2d(1,11))
      CALL SUM_K_V( NNstat, kmax, m_rho_w_z1_2, wrk2d(1,3), &
           wrk2d(1,11))

      DO j = 1,NNstat
         MS_RUSS(j) = MS_RUSS(j) + wrk2d(j,1)
         MS_RVSS(j) = MS_RVSS(j) + wrk2d(j,2)
         MS_RWSS(j) = MS_RWSS(j) + wrk2d(j,3)
      ENDDO

#undef m_rho_u_z1_2
#undef m_rho_v_z1_2
#undef m_rho_w_z1_2

#undef m_rho_z1

#define m_z1_z1 tmp1
#define m_z1_u  tmp2
#define m_z1_v  tmp3
#define m_z1_w  tmp4

      DO j = 1,NNstat*kmax
         m_z1_z1(j) = m_z1(j)*m_z1(j)
         m_z1_u(j) = m_z1(j)*m_u(j)
         m_z1_v(j) = m_z1(j)*m_v(j)
         m_z1_w(j) = m_z1(j)*m_w(j)
      ENDDO
      
      CALL SUM_K_V( NNstat, kmax, m_z1_z1, wrk2d(1,1), &
           wrk2d(1,11) )
      CALL SUM_K_V( NNstat, kmax, m_z1_u, wrk2d(1,2), &
           wrk2d(1,11) )
      CALL SUM_K_V( NNstat, kmax, m_z1_v, wrk2d(1,3), &
           wrk2d(1,11) )
      CALL SUM_K_V( NNstat, kmax, m_z1_w, wrk2d(1,4), &
           wrk2d(1,11) )

      DO j = 1,NNstat
         MS_S2(j) = MS_S2(j) + wrk2d(j,1)
         MS_SU(j) = MS_SU(j) + wrk2d(j,2)
         MS_SV(j) = MS_SV(j) + wrk2d(j,3)
         MS_SW(j) = MS_SW(j) + wrk2d(j,4)
      ENDDO

#undef m_z1_z1 
#undef m_z1_u  
#undef m_z1_v  
#undef m_z1_w

! # Intermittency

#define m_z1_gamma tmp1

#ifdef SINGLE_PREC
      inter1 = 0.02e0
#else
      inter1 = 0.02d0
#endif
      inter2 = C_1_R-inter1

      DO j = 1,NNstat*kmax
         IF ( inter1 .LE. m_z1(j) .AND. m_z1(j) .LE. inter2 ) THEN
            m_z1_gamma(j) = C_1_R
         ELSE
            m_z1_gamma(j) = C_0_R
         ENDIF
      ENDDO

      CALL SUM_K_V( NNstat, kmax, m_z1_gamma, wrk2d(1,1), &
           wrk2d(1,11) )

      DO j = 1,NNstat
         MS_GAMMA(j) = MS_GAMMA(j) + wrk2d(j,1)
      ENDDO

#undef m_z1_gamma

! # Skewness and flatness

#define m_z1_3 tmp1
#define m_z1_4 tmp2

      DO j = 1,NNstat*kmax
         m_z1_3(j) = m_z1(j)*m_z1(j)*m_z1(j)
         m_z1_4(j) = m_z1_3(j)*m_z1(j)
      ENDDO
      CALL SUM_K_V( NNstat, kmax, m_z1_3, wrk2d(1,1), &
           wrk2d(1,11) )
      CALL SUM_K_V( NNstat, kmax, m_z1_4, wrk2d(1,2), &
           wrk2d(1,11) )
      DO j = 1,NNstat
         MS_S3(j) = MS_S3(j) + wrk2d(j,1)
         MS_S4(j) = MS_S4(j) + wrk2d(j,2)
      ENDDO

#undef m_z1_3
#undef m_z1_4

      RETURN
      END

