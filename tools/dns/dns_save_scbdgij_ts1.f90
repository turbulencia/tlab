#include "types.h"
#include "dns_const.h"
#include "avgij_map.h"

SUBROUTINE DNS_SAVE_SCBDGIJ_TS1(NNstat, m_z1, u, v, w, vis, z1, tmp1, tmp2, tmp3, tmp4, tmp5, &
     tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, mean1d_sc, wrk2d, wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : nstatavg, statavg
  USE DNS_GLOBAL, ONLY : isize_wrk2d
  USE DNS_GLOBAL, ONLY : itransport, visc

  IMPLICIT NONE

  TINTEGER NNstat
  TREAL m_z1(*)
  TREAL tmp1(*)
  TREAL tmp2(*)
  TREAL tmp3(*)
  TREAL tmp4(*)
  TREAL tmp5(*)
  TREAL tmp6(*)
  TREAL tmp7(*)
  TREAL tmp8(*)
  TREAL tmp9(*)
  TREAL tmp10(*)
  TREAL tmp11(*)

  TREAL u(*)
  TREAL v(*)
  TREAL w(*)
  TREAL vis(*)
  TREAL z1(imax,jmax,kmax)
  TREAL mean1d_sc(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)
  TREAL wrk3d(*)
  TINTEGER j, bcs(2,1)
  TREAL c23, c43, aux1, aux2

  c23 = C_2_R/C_3_R
  c43 = C_4_R/C_3_R

  bcs = 0

  CALL REDUCE( imax,jmax,kmax, z1, nstatavg, statavg, m_z1)

#define m_vis  tmp5
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL REDUCE(imax,jmax,kmax, vis, nstatavg, statavg, m_vis)
  ELSE
     DO j = 1,NNstat*kmax
        m_vis(j) = C_1_R
     ENDDO
  ENDIF

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1, tmp6, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1, tmp7, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1, tmp8, wrk3d, wrk2d,wrk3d)

#define m_z1_x tmp9
#define m_z1_y tmp10
#define m_z1_z tmp11

  CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_z1_x)
  CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_z1_y)
  CALL REDUCE(imax,jmax,kmax, tmp8, nstatavg, statavg, m_z1_z)

! Cross terms xy
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp6, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp7, wrk3d, wrk2d,wrk3d)

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, kmax*imax*jmax
        tmp8(j) = vis(j)*(tmp6(j)+tmp7(j))
     ENDDO
  ELSE
     DO j=1, kmax*imax*jmax
        tmp8(j) = tmp6(j)+tmp7(j)
     ENDDO
  ENDIF

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp8, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp8, tmp2, wrk3d, wrk2d,wrk3d)

#define m_v_x tmp8
#define m_u_y tmp6

  CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_v_x)
  CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_u_y)

#define m_tau_xy_y tmp3
#define m_tau_xy_x tmp4

  CALL REDUCE(imax,jmax,kmax, tmp1, nstatavg, statavg, m_tau_xy_y)
  CALL REDUCE(imax,jmax,kmax, tmp2, nstatavg, statavg, m_tau_xy_x)

#define m_z1_tau_xy_y tmp1
#define m_z1_tau_xy_x tmp2
#define m_fx_u_y      tmp7
#define m_fx_v_x      wrk3d
#define m_fy_u_y      tmp6
#define m_fy_v_x      tmp8

  DO j = 1,NNstat*kmax
     aux1 = m_u_y(j)
     aux2 = m_v_x(j)
     m_fx_u_y(j) = m_vis(j)*m_z1_x(j)*aux1
     m_fx_v_x(j) = m_vis(j)*m_z1_x(j)*aux2
     m_fy_u_y(j) = m_vis(j)*m_z1_y(j)*aux1
     m_fy_v_x(j) = m_vis(j)*m_z1_y(j)*aux2
     m_z1_tau_xy_y(j) = m_tau_xy_y(j)*m_z1(j)
     m_z1_tau_xy_x(j) = m_tau_xy_x(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_fx_u_y, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fx_v_x, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fy_u_y, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fy_v_x, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_tau_xy_y, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_xy_y, wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_xy_x, wrk2d(1,7), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_tau_xy_x, wrk2d(1,8), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkVk(j) = MS_FkVk(j) + wrk2d(j,2)
     MS_FkUk(j) = MS_FkUk(j) + wrk2d(j,3)

     MS_TAUxkSk(j) = MS_TAUxkSk(j)    + visc*(wrk2d(j,3) + wrk2d(j,4))
     MS_TAUykSk(j) = MS_TAUykSk(j)    + visc*(wrk2d(j,1) + wrk2d(j,2))

     MS_TAUxyy(j) = MS_TAUxyy(j) + visc*wrk2d(j,5)
     MS_STAUxkk(j) = MS_STAUxkk(j) + visc*wrk2d(j,6)
     MS_TAUxyx(j) = MS_TAUxyx(j) + visc*wrk2d(j,8)
     MS_STAUykk(j) = MS_STAUykk(j) + visc*wrk2d(j,7)
  ENDDO

#undef m_fx_u_y 
#undef m_fx_v_x 
#undef m_fy_u_y 
#undef m_fy_v_x 
#undef m_v_x 
#undef m_u_y 
#undef m_tau_xy_y
#undef m_tau_xy_x
#undef m_z1_tau_xy_y
#undef m_z1_tau_xy_x

! Cross terms xz

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tmp6, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, tmp7, wrk3d, wrk2d,wrk3d)

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, kmax*imax*jmax
        tmp8(j) = vis(j)*(tmp6(j)+tmp7(j))
     ENDDO
  ELSE
     DO j=1, kmax*imax*jmax
        tmp8(j) = tmp6(j)+tmp7(j)
     ENDDO
  ENDIF

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp8, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp8, tmp2, wrk3d, wrk2d,wrk3d)

#define m_w_x tmp8
#define m_u_z tmp6

  CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_w_x)
  CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_u_z)

#define m_tau_xz_z tmp3
#define m_tau_xz_x tmp4

  CALL REDUCE(imax,jmax,kmax, tmp1, nstatavg, statavg, m_tau_xz_z)
  CALL REDUCE(imax,jmax,kmax, tmp2, nstatavg, statavg, m_tau_xz_x)

#define m_z1_tau_xz_z tmp1
#define m_z1_tau_xz_x tmp2
#define m_fx_u_z      tmp7
#define m_fx_w_x      wrk3d
#define m_fz_u_z      tmp8
#define m_fz_w_x      tmp6

  DO j = 1,NNstat*kmax
     aux1 = m_u_z(j)
     aux2 = m_w_x(j)
     m_fx_u_z(j) = m_vis(j)*m_z1_x(j)*aux1
     m_fx_w_x(j) = m_vis(j)*m_z1_x(j)*aux2
     m_fz_u_z(j) = m_vis(j)*m_z1_z(j)*aux1
     m_fz_w_x(j) = m_vis(j)*m_z1_z(j)*aux2
     m_z1_tau_xz_z(j) = m_tau_xz_z(j)*m_z1(j)
     m_z1_tau_xz_x(j) = m_tau_xz_x(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_fx_u_z, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fx_w_x, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_u_z, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_w_x, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_tau_xz_x, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_xz_z, wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_xz_x, wrk2d(1,7), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkWk(j) = MS_FkWk(j) + wrk2d(j,2)
     MS_FkUk(j) = MS_FkUk(j) + wrk2d(j,3)

     MS_TAUxkSk(j) = MS_TAUxkSk(j)    + visc*(wrk2d(j,3) + wrk2d(j,4))
     MS_TAUzkSk(j) = MS_TAUzkSk(j)    + visc*(wrk2d(j,1) + wrk2d(j,2))

     MS_STAUxkk(j) = MS_STAUxkk(j) + visc*wrk2d(j,6)
     MS_TAUxzx(j) = MS_TAUxzx(j) + visc*wrk2d(j,5)
     MS_STAUzkk(j) = MS_STAUzkk(j) + visc*wrk2d(j,7)
  ENDDO

#undef m_fx_u_z
#undef m_fx_w_x 
#undef m_fz_u_z
#undef m_fz_w_x 
#undef m_w_x 
#undef m_u_z
#undef m_tau_xz_z
#undef m_tau_xz_x
#undef m_z1_tau_xz_z
#undef m_z1_tau_xz_x

! Cross terms yz

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp6, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, tmp7, wrk3d, wrk2d,wrk3d)

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, kmax*imax*jmax
        tmp8(j) = vis(j)*(tmp6(j)+tmp7(j))
     ENDDO
  ELSE
     DO j=1, kmax*imax*jmax
        tmp8(j) = tmp6(j)+tmp7(j)
     ENDDO
  ENDIF

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp8, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp8, tmp2, wrk3d, wrk2d,wrk3d)

#define m_w_y tmp8
#define m_v_z tmp6

  CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_w_y)
  CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_v_z)

#define m_tau_yz_z tmp3
#define m_tau_yz_y tmp4

  CALL REDUCE(imax,jmax,kmax, tmp1, nstatavg, statavg, m_tau_yz_z)
  CALL REDUCE(imax,jmax,kmax, tmp2, nstatavg, statavg, m_tau_yz_y)

#define m_z1_tau_yz_z tmp1
#define m_z1_tau_yz_y tmp2
#define m_fy_v_z      tmp7
#define m_fy_w_y      wrk3d
#define m_fz_v_z      tmp8
#define m_fz_w_y      tmp6

  DO j = 1,NNstat*kmax
     aux1 = m_v_z(j)
     aux2 = m_w_y(j)
     m_fy_v_z(j) = m_vis(j)*m_z1_y(j)*aux1
     m_fy_w_y(j) = m_vis(j)*m_z1_y(j)*aux2
     m_fz_v_z(j) = m_vis(j)*m_z1_z(j)*aux1
     m_fz_w_y(j) = m_vis(j)*m_z1_z(j)*aux2
     m_z1_tau_yz_z(j) = m_tau_yz_z(j)*m_z1(j)
     m_z1_tau_yz_y(j) = m_tau_yz_y(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_fy_v_z, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fy_w_y, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_v_z, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_w_y, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_tau_yz_y, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_yz_z, wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_yz_y, wrk2d(1,7), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkWk(j) = MS_FkWk(j) + wrk2d(j,2)
     MS_FkVk(j) = MS_FkVk(j) + wrk2d(j,3)

     MS_TAUykSk(j) = MS_TAUykSk(j)    + visc*(wrk2d(j,3) + wrk2d(j,4))
     MS_TAUzkSk(j) = MS_TAUzkSk(j)    + visc*(wrk2d(j,1) + wrk2d(j,2))

     MS_TAUyzy(j) = MS_TAUyzy(j) + visc*wrk2d(j,5)
     MS_STAUykk(j) = MS_STAUykk(j) + visc*wrk2d(j,6)
     MS_STAUzkk(j) = MS_STAUzkk(j) + visc*wrk2d(j,7)
  ENDDO

#undef m_fy_v_z
#undef m_fy_w_y 
#undef m_fz_v_z
#undef m_fz_w_y 
#undef m_w_y 
#undef m_v_z
#undef m_tau_yz_z
#undef m_tau_yz_y
#undef m_z1_tau_yz_z
#undef m_z1_tau_yz_y

! Cross terms xx

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)

#define m_u_x tmp8

  CALL REDUCE(imax,jmax,kmax, tmp1, nstatavg, statavg, m_u_x)

#define m_fx_u_x tmp7
#define m_fy_u_x wrk3d
#define m_fz_u_x tmp6

  DO j = 1,NNstat*kmax
     m_fx_u_x(j) = m_vis(j)*m_z1_x(j)*m_u_x(j)
     m_fy_u_x(j) = m_vis(j)*m_z1_y(j)*m_u_x(j)
     m_fz_u_x(j) = m_vis(j)*m_z1_z(j)*m_u_x(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_fx_u_x, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fy_u_x, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_u_x, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkUk(j) = MS_FkUk(j) + wrk2d(j,1)
     MS_TAUxkSk(j) = MS_TAUxkSk(j)    + c43*wrk2d(j,1)*visc
     MS_TAUykSk(j) = MS_TAUykSk(j)    - c23*wrk2d(j,2)*visc
     MS_TAUzkSk(j) = MS_TAUzkSk(j)    - c23*wrk2d(j,3)*visc
  ENDDO

#undef m_fx_u_x
#undef m_fy_u_x 
#undef m_fz_u_x
#undef m_u_x

! Cross terms yy

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)

#define m_v_y tmp8
  CALL REDUCE(imax,jmax,kmax, tmp2, nstatavg, statavg, m_v_y)

#define m_fx_v_y tmp7
#define m_fy_v_y wrk3d
#define m_fz_v_y tmp6

  DO j = 1,NNstat*kmax
     m_fx_v_y(j) = m_vis(j)*m_z1_x(j)*m_v_y(j)
     m_fy_v_y(j) = m_vis(j)*m_z1_y(j)*m_v_y(j)
     m_fz_v_y(j) = m_vis(j)*m_z1_z(j)*m_v_y(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_fx_v_y, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fy_v_y, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_v_y, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkVk(j) = MS_FkVk(j) + wrk2d(j,2)
     MS_TAUykSk(j) = MS_TAUykSk(j)    + c43*wrk2d(j,2)*visc
     MS_TAUxkSk(j) = MS_TAUxkSk(j)    - c23*wrk2d(j,1)*visc
     MS_TAUzkSk(j) = MS_TAUzkSk(j)    - c23*wrk2d(j,3)*visc
  ENDDO

#undef m_fx_v_y
#undef m_fy_v_y 
#undef m_fz_v_y
#undef m_v_y

! Cross terms zz

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)

#define m_w_z tmp8
  CALL REDUCE(imax,jmax,kmax, tmp3, nstatavg, statavg, m_w_z)

#define m_fx_w_z tmp7
#define m_fy_w_z wrk3d
#define m_fz_w_z tmp6

  DO j = 1,NNstat*kmax
     m_fx_w_z(j) = m_vis(j)*m_z1_x(j)*m_w_z(j)
     m_fy_w_z(j) = m_vis(j)*m_z1_y(j)*m_w_z(j)
     m_fz_w_z(j) = m_vis(j)*m_z1_z(j)*m_w_z(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_fx_w_z, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fy_w_z, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_fz_w_z, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_FkWk(j) = MS_FkWk(j) + wrk2d(j,3)
     MS_TAUzkSk(j) = MS_TAUzkSk(j)    + c43*wrk2d(j,3)*visc
     MS_TAUxkSk(j) = MS_TAUxkSk(j)    - c23*wrk2d(j,1)*visc
     MS_TAUykSk(j) = MS_TAUykSk(j)    - c23*wrk2d(j,2)*visc
  ENDDO

#undef m_fx_w_z
#undef m_fy_w_z 
#undef m_fz_w_z
#undef m_w_z

#undef m_vis  

! Compute Cross terms

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j=1, imax*jmax*kmax
        tmp4(j) = vis(j)*(c43*tmp1(j)-c23*(tmp2(j)+tmp3(j)))
        tmp5(j) = vis(j)*(c43*tmp2(j)-c23*(tmp1(j)+tmp3(j)))
        tmp6(j) = vis(j)*(c43*tmp3(j)-c23*(tmp1(j)+tmp2(j)))
     ENDDO
  ELSE
     DO j=1, imax*jmax*kmax
        tmp4(j) = c43*tmp1(j)-c23*(tmp2(j)+tmp3(j))
        tmp5(j) = c43*tmp2(j)-c23*(tmp1(j)+tmp3(j))
        tmp6(j) = c43*tmp3(j)-c23*(tmp1(j)+tmp2(j))
     ENDDO
  ENDIF

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp6, tmp3, wrk3d, wrk2d,wrk3d)

#define m_tau_xx_x tmp4
#define m_tau_yy_y tmp5
#define m_tau_zz_z tmp6

  CALL REDUCE(imax,jmax,kmax, tmp1, nstatavg, statavg, m_tau_xx_x)
  CALL REDUCE(imax,jmax,kmax, tmp2, nstatavg, statavg, m_tau_yy_y)
  CALL REDUCE(imax,jmax,kmax, tmp3, nstatavg, statavg, m_tau_zz_z)

#define m_z1_tau_xx_x tmp1
#define m_z1_tau_yy_y tmp2
#define m_z1_tau_zz_z tmp3

  DO j = 1,NNstat*kmax
     m_z1_tau_xx_x(j) = m_tau_xx_x(j)*m_z1(j)
     m_z1_tau_yy_y(j) = m_tau_yy_y(j)*m_z1(j)
     m_z1_tau_zz_z(j) = m_tau_zz_z(j)*m_z1(j)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, m_tau_xx_x, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_tau_yy_y, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_xx_x, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_yy_y, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, m_z1_tau_zz_z, wrk2d(1,5), wrk2d(1,11) )

  DO j = 1,NNstat
     MS_TAUxxx(j) = MS_TAUxxx(j) + wrk2d(j,1)*visc
     MS_TAUyyy(j) = MS_TAUyyy(j) + wrk2d(j,2)*visc
     MS_STAUxkk(j) = MS_STAUxkk(j) + wrk2d(j,3)*visc
     MS_STAUykk(j) = MS_STAUykk(j) + wrk2d(j,4)*visc
     MS_STAUzkk(j) = MS_STAUzkk(j) + wrk2d(j,5)*visc
  ENDDO

#undef m_tau_xx_x 
#undef m_tau_yy_y 
#undef m_tau_zz_z 
#undef m_z1_tau_xx_x 
#undef m_z1_tau_yy_y 
#undef m_z1_tau_zz_z 

#undef m_z1_x 
#undef m_z1_y 
#undef m_z1_z 

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ_TS1

