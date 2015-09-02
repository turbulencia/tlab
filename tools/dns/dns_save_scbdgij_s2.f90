#include "types.h"
#include "dns_const.h"
#include "avgij_map.h"

SUBROUTINE DNS_SAVE_SCBDGIJ_S2(NNstat, dx, dy, dz, &
     m_rho, m_u, m_v, m_w, m_z1, p, vis, z1, &
     tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, &
     mean1d_sc, wrk1d, wrk2d, wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : imode_fdm, i1bc,j1bc,k1bc
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d, itransport

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL dx(*)
  TREAL dy(*)
  TREAL dz(*)
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
  TREAL tmp7(*)

  TREAL p(*)
  TREAL vis(*)
  TREAL z1(imax,jmax,kmax)
  TREAL mean1d_sc(nstatavg,jmax,*)
  TREAL wrk1d(*)
  TREAL wrk2d(isize_wrk2d,*)
  TREAL wrk3d(*)

  TINTEGER j, i0

  i0 = 0

  CALL REDUCE( imax, jmax, kmax, z1, nstatavg, statavg, m_z1)

  CALL PARTIAL_X( imode_fdm, imax, jmax, kmax, i1bc, &
       dx, z1, tmp2, i0, i0, wrk1d, wrk2d, wrk3d )
  CALL PARTIAL_Y( imode_fdm, imax, jmax, kmax, j1bc,&
       dy, z1, tmp3, i0, i0, wrk1d, wrk2d, wrk3d )
  CALL PARTIAL_Z( imode_fdm, imax, jmax, kmax, k1bc,&
       dz, z1, tmp4, i0, i0, wrk1d, wrk2d, wrk3d )

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, kmax*imax*jmax
        tmp5(j) = vis(j)*tmp2(j)
     ENDDO
  ELSE
     DO j=1, kmax*imax*jmax
        tmp5(j) = tmp2(j)
     ENDDO
  ENDIF

  CALL PARTIAL_X( imode_fdm, imax, jmax, kmax, i1bc, &
       dx, tmp5, tmp6, i0, i0, wrk1d, wrk2d, wrk3d )

#define m_fxx tmp7

  CALL REDUCE(imax, jmax, kmax, tmp6, nstatavg, statavg, m_fxx)

  CALL SUM_K_V( NNstat, kmax, m_fxx, wrk2d(1,1), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_Fxx(j) = MS_Fxx(j) + wrk2d(j,1)
  ENDDO

#define m_u_fxx tmp1
#define m_v_fxx wrk3d
#define m_w_fxx tmp5
#define m_z1_fxx m_fxx

  DO j = 1,NNstat*kmax
     m_u_fxx(j) = m_fxx(j)*m_u(j)
     m_v_fxx(j) = m_fxx(j)*m_v(j)
     m_w_fxx(j) = m_fxx(j)*m_w(j)
     m_z1_fxx(j) = m_fxx(j)*m_z1(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_z1_fxx, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_u_fxx, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_v_fxx, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_w_fxx, wrk2d(1,4), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkdkS(j) = MS_FkdkS(j) + wrk2d(j,1)
     MS_FkdkU(j) = MS_FkdkU(j) + wrk2d(j,2)
     MS_FkdkV(j) = MS_FkdkV(j) + wrk2d(j,3)
     MS_FkdkW(j) = MS_FkdkW(j) + wrk2d(j,4)
  ENDDO

#undef m_u_fxx 
#undef m_v_fxx 
#undef m_w_fxx 
#undef m_fxx
#undef m_z1_fxx

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, kmax*imax*jmax
        tmp5(j) = vis(j)*tmp3(j)
     ENDDO
  ELSE
     DO j=1, kmax*imax*jmax
        tmp5(j) = tmp3(j)
     ENDDO
  ENDIF
  CALL PARTIAL_Y( imode_fdm, imax, jmax, kmax, j1bc, &
       dy, tmp5, tmp6, i0, i0, wrk1d, wrk2d, wrk3d )

#define m_fyy tmp7

  CALL REDUCE(imax, jmax, kmax, tmp6, nstatavg, statavg, m_fyy)

  CALL SUM_K_V( NNstat, kmax, m_fyy, wrk2d(1,1), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_Fyy(j) = MS_Fyy(j) + wrk2d(j,1)
  ENDDO

#define m_u_fyy tmp1
#define m_v_fyy wrk3d
#define m_w_fyy tmp5
#define m_z1_fyy m_fyy

  DO j = 1,NNstat*kmax
     m_u_fyy(j) = m_fyy(j)*m_u(j)
     m_v_fyy(j) = m_fyy(j)*m_v(j)
     m_w_fyy(j) = m_fyy(j)*m_w(j)
     m_z1_fyy(j) = m_fyy(j)*m_z1(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_z1_fyy, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_u_fyy, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_v_fyy, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_w_fyy, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkdkS(j) = MS_FkdkS(j) + wrk2d(j,1)
     MS_FkdkU(j) = MS_FkdkU(j) + wrk2d(j,2)
     MS_FkdkV(j) = MS_FkdkV(j) + wrk2d(j,3)
     MS_FkdkW(j) = MS_FkdkW(j) + wrk2d(j,4)
  ENDDO

#undef m_u_fyy 
#undef m_v_fyy 
#undef m_w_fyy 
#undef m_fyy
#undef m_z1_fyy

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, kmax*imax*jmax
        tmp5(j) = vis(j)*tmp4(j)
     ENDDO
  ELSE
     DO j=1, kmax*imax*jmax
        tmp5(j) = tmp4(j)
     ENDDO
  ENDIF

  CALL PARTIAL_Z( imode_fdm, imax, jmax, kmax, k1bc, &
       dz, tmp5, tmp6, i0, i0, wrk1d, wrk2d, wrk3d )

#define m_fzz tmp7

  CALL REDUCE(imax, jmax, kmax, tmp6, nstatavg, statavg, m_fzz)

#define m_u_fzz tmp1
#define m_v_fzz wrk3d
#define m_w_fzz tmp5
#define m_z1_fzz m_fzz

  DO j = 1,NNstat*kmax
     m_u_fzz(j) = m_fzz(j)*m_u(j)
     m_v_fzz(j) = m_fzz(j)*m_v(j)
     m_w_fzz(j) = m_fzz(j)*m_w(j)
     m_z1_fzz(j) = m_fzz(j)*m_z1(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_z1_fzz, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_u_fzz, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_v_fzz, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_w_fzz, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkdkS(j) = MS_FkdkS(j) + wrk2d(j,1)
     MS_FkdkU(j) = MS_FkdkU(j) + wrk2d(j,2)
     MS_FkdkV(j) = MS_FkdkV(j) + wrk2d(j,3)
     MS_FkdkW(j) = MS_FkdkW(j) + wrk2d(j,4)
  ENDDO

#undef m_u_fzz 
#undef m_v_fzz 
#undef m_w_fzz 
#undef m_fzz
#undef m_z1_fzz

#define m_z1_x tmp5
#define m_z1_y tmp6
#define m_z1_z tmp7

  CALL REDUCE(imax, jmax, kmax, tmp2, nstatavg, statavg, m_z1_x)
  CALL REDUCE(imax, jmax, kmax, tmp3, nstatavg, statavg, m_z1_y)
  CALL REDUCE(imax, jmax, kmax, tmp4, nstatavg, statavg, m_z1_z)

#define m_vis  tmp1
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL REDUCE(imax, jmax, kmax, vis, nstatavg, statavg, m_vis)      
  ELSE
     DO j = 1,NNstat*kmax
        m_vis(j) = C_1_R
     ENDDO
  ENDIF

  CALL SUM_K_V( NNstat, kmax, m_z1_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_z1_y, wrk2d(1,2), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_Sx(j) = MS_Sx(j) + wrk2d(j,1)
     MS_Sy(j) = MS_Sy(j) + wrk2d(j,2)
  ENDDO

#define m_vis_z1_x tmp2
#define m_vis_z1_y tmp3
#define m_vis_z1_z tmp4

  DO j = 1,NNstat*kmax
     m_vis_z1_x(j) = m_vis(j)*m_z1_x(j)
     m_vis_z1_y(j) = m_vis(j)*m_z1_y(j)
     m_vis_z1_z(j) = m_vis(j)*m_z1_z(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_vis_z1_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_vis_z1_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_vis_z1_z, wrk2d(1,3), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_Fx(j) = MS_Fx(j) + wrk2d(j,1)
     MS_Fy(j) = MS_Fy(j) + wrk2d(j,2)
     MS_Fz(j) = MS_Fz(j) + wrk2d(j,3)
  ENDDO

! here m_vis_z1_z contains the total dissipation
#define m_eps m_vis_z1_z
  DO j = 1,NNstat*kmax
     m_eps(j) = m_vis(j)*(m_z1_x(j)**2&
          +m_z1_y(j)**2+m_z1_z(j)**2)
     m_vis_z1_x(j) = m_vis_z1_x(j)*m_z1(j)
     m_vis_z1_y(j) = m_vis_z1_y(j)*m_z1(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_vis_z1_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_vis_z1_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_eps, wrk2d(1,3), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_SFx(j) = MS_SFx(j) + wrk2d(j,1)
     MS_SFy(j) = MS_SFy(j) + wrk2d(j,2)
     MS_SEPS(j) = MS_SEPS(j) + wrk2d(j,3)
  ENDDO

#undef m_vis_z1_x 
#undef m_vis_z1_y 
#undef m_vis_z1_z 
#undef m_eps
#undef m_vis

#define m_p   tmp1
  CALL REDUCE( imax, jmax, kmax, p,  nstatavg, statavg, m_p)

#define m_p_z1_x tmp3
#define m_p_z1_y wrk3d
#define m_p_z1_z tmp2

  DO j = 1,NNstat*kmax
     m_p_z1_x(j) = m_p(j)*m_z1_x(j)
     m_p_z1_y(j) = m_p(j)*m_z1_y(j)
     m_p_z1_z(j) = m_p(j)*m_z1_z(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_p_z1_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_p_z1_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_p_z1_z, wrk2d(1,3), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_PSx(j) = MS_PSx(j) + wrk2d(j,1)
     MS_PSy(j) = MS_PSy(j) + wrk2d(j,2)
     MS_PSz(j) = MS_PSz(j) + wrk2d(j,3)
  ENDDO

#undef m_p_z1_x
#undef m_p_z1_y
#undef m_p_z1_z

#define m_rho_z1_x    tmp1
#define m_rho_z1_y    tmp2
#define m_rho_z1_z1_x tmp3
#define m_rho_z1_z1_y tmp4

  DO j = 1,NNstat*kmax
     m_rho_z1_x(j) = m_rho(j)*m_z1_x(j)
     m_rho_z1_y(j) = m_rho(j)*m_z1_y(j)
     m_rho_z1_z1_x(j) = m_rho_z1_x(j)*m_z1(j)
     m_rho_z1_z1_y(j) = m_rho_z1_y(j)*m_z1(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_rho_z1_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_z1_x, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_z1_y, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSx(j) = MS_RSx(j) + wrk2d(j,1)
     MS_RSy(j) = MS_RSy(j) + wrk2d(j,2)
     MS_RSSx(j) = MS_RSSx(j) + C_2_R*wrk2d(j,3)
     MS_RSSy(j) = MS_RSSy(j) + C_2_R*wrk2d(j,4)
  ENDDO

#undef m_rho_z1_x    
#undef m_rho_z1_y    
#undef m_rho_z1_z1_x 
#undef m_rho_z1_z1_y 

#define m_rho_z1_x_u    tmp1
#define m_rho_z1_y_v    tmp2
#define m_rho_z1_z1_x_u tmp3
#define m_rho_z1_z1_y_v tmp4

  DO j = 1,NNstat*kmax
     m_rho_z1_x_u(j) = m_rho(j)*m_z1_x(j)*m_u(j)
     m_rho_z1_y_v(j) = m_rho(j)*m_z1_y(j)*m_v(j)
     m_rho_z1_z1_x_u(j) = m_rho_z1_x_u(j)*m_z1(j)
     m_rho_z1_z1_y_v(j) = m_rho_z1_y_v(j)*m_z1(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_rho_z1_x_u, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_y_v, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_z1_x_u, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_z1_y_v, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSUx(j) = MS_RSUx(j) + wrk2d(j,1)
     MS_RSVy(j) = MS_RSVy(j) + wrk2d(j,2)
     MS_RSSUx(j) = MS_RSSUx(j) + C_2_R*wrk2d(j,3)
     MS_RSSVy(j) = MS_RSSVy(j) + C_2_R*wrk2d(j,4)
  ENDDO

#undef m_rho_z1_x_u
#undef m_rho_z1_y_v
#undef m_rho_z1_z1_x_u 
#undef m_rho_z1_z1_y_v 

#define m_rho_z1_x_v tmp1
#define m_rho_z1_y_u tmp2
#define m_rho_z1_x_w tmp3
#define m_rho_z1_y_w tmp4

  DO j = 1,NNstat*kmax
     m_rho_z1_x_v(j) = m_rho(j)*m_z1_x(j)*m_v(j)
     m_rho_z1_y_u(j) = m_rho(j)*m_z1_y(j)*m_u(j)
     m_rho_z1_x_w(j) = m_rho(j)*m_z1_x(j)*m_w(j)
     m_rho_z1_y_w(j) = m_rho(j)*m_z1_y(j)*m_w(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_rho_z1_x_v, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_y_u, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_x_w, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_z1_y_w, wrk2d(1,4), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RSVx(j) = MS_RSVx(j) + wrk2d(j,1)
     MS_RSUy(j) = MS_RSUy(j) + wrk2d(j,2)
     MS_RSWx(j) = MS_RSWx(j) + wrk2d(j,3)
     MS_RSWy(j) = MS_RSWy(j) + wrk2d(j,4)
  ENDDO

#undef m_rho_z1_x_v
#undef m_rho_z1_y_u
#undef m_rho_z1_x_w 
#undef m_rho_z1_y_w


#define m_rho_u_u_z1_x tmp1
#define m_rho_v_v_z1_y tmp2
#define m_rho_u_v_z1_x tmp3
#define m_rho_u_v_z1_y tmp4
#define m_rho_u_w_z1_x wrk3d
#define m_rho_v_w_z1_y tmp7

  DO j = 1,NNstat*kmax
     m_rho_u_u_z1_x(j) = m_rho(j)*m_z1_x(j)*m_u(j)**2
     m_rho_v_v_z1_y(j) = m_rho(j)*m_z1_y(j)*m_v(j)**2
     m_rho_u_v_z1_x(j) = m_rho(j)*m_z1_x(j)*m_u(j)*m_v(j)
     m_rho_u_v_z1_y(j) = m_rho(j)*m_z1_y(j)*m_u(j)*m_v(j)
     m_rho_u_w_z1_x(j) = m_rho(j)*m_z1_x(j)*m_u(j)*m_w(j)
     m_rho_v_w_z1_y(j) = m_rho(j)*m_z1_y(j)*m_v(j)*m_w(j)
  ENDDO

  CALL SUM_K_V( NNstat, kmax, m_rho_u_u_z1_x, wrk2d(1,1), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_v_v_z1_y, wrk2d(1,2), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_u_v_z1_x, wrk2d(1,3), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_u_v_z1_y, wrk2d(1,4), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_u_w_z1_x, wrk2d(1,5), &
       wrk2d(1,11) )
  CALL SUM_K_V( NNstat, kmax, m_rho_v_w_z1_y, wrk2d(1,6), &
       wrk2d(1,11) )

  DO j = 1,NNstat
     MS_RUUSx(j) = MS_RUUSx(j) + wrk2d(j,1)
     MS_RVVSy(j) = MS_RVVSy(j) + wrk2d(j,2)
     MS_RUVSx(j) = MS_RUVSx(j) + wrk2d(j,3)
     MS_RUVSy(j) = MS_RUVSy(j) + wrk2d(j,4)
     MS_RUWSx(j) = MS_RUWSx(j) + wrk2d(j,5)
     MS_RVWSy(j) = MS_RVWSy(j) + wrk2d(j,6)
  ENDDO

#undef m_rho_u_u_z1_x
#undef m_rho_v_v_z1_y
#undef m_rho_u_v_z1_x
#undef m_rho_u_v_z1_y
#undef m_rho_u_w_z1_x
#undef m_rho_v_w_z1_y

#undef m_z1_x
#undef m_z1_y
#undef m_z1_z

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_S2

