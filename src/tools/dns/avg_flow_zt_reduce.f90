#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "avgij_map.h"

! #####################################################
! # Preparing data for statistics.
! # Clean-up means the temporary arrays can be
! # over-written.
! #
! # Lay-out:
! #    Single terms
! #    Double-cross terms
! #    Triple-cross terms
! #    Scalar terms
! #    Derivatives in u
! #    Derivatives in v
! #    Derivatives in w
! #    Shear Stress tensor and its derivatives
! #    Derivatives in rho
! #    Derivatives in p
! #    Derivatives in T
! #
! # 09/25/00 Juan Pedro Mellado
! #####################################################
SUBROUTINE AVG_FLOW_ZT_REDUCE(q, hq,txc, mean1d, wrk2d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
#endif
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_wrk2d, imode_eqns
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : itransport, visc
  USE TLAB_VARS, ONLY : nstatavg, statavg, nstatavg_points
  USE TLAB_PROCS
  USE AVGS, ONLY: SUM1V1D_V
  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN   ), TARGET :: q
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(INOUT), TARGET :: txc, hq
  TREAL mean1d(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)
  TREAL wrk3d(imax,jmax,kmax)

  TINTEGER j, bcs(2,1)
  TINTEGER NNstat
  TREAL c2, c23, cs2

  ! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u, v, w, rho, p, T, vis
  TREAL, DIMENSION(:,:,:), POINTER :: xc, yc, zc, vc, wc, uc, tc, sc, rc, qc, oc, pc

  ! ###################################################################
#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'ENTERING AVG_FLOW_ZT_REDUCE' )
#endif

  IF ( imax .LT. nstatavg ) THEN
    CALL TLAB_WRITE_ASCII(efile, 'AVG_FLOW_ZT_REDUCE. Not enough space in available arrays.')
    CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  nstatavg_points = nstatavg_points + g(3)%size

  bcs = 0
  c2 = C_2_R
  c23 = C_2_R/C_3_R
  cs2 = C_14_R*C_1EM2_R
  NNstat = nstatavg*jmax

  ! Define pointers
  u   => q(:,:,:,1)
  v   => q(:,:,:,2)
  w   => q(:,:,:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     rho => q(:,:,:,5)
     p   => q(:,:,:,6)
     T   => q(:,:,:,7)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
  ENDIF

  xc => hq(:,:,:,1)
  yc => hq(:,:,:,2)
  zc => hq(:,:,:,3)
  vc => txc(:,:,:,1)
  wc => txc(:,:,:,2)
  uc => txc(:,:,:,3)
  tc => txc(:,:,:,4)
  sc => txc(:,:,:,5)
  rc => txc(:,:,:,6)
  qc => txc(:,:,:,7)
  oc => txc(:,:,:,8)
  pc => txc(:,:,:,9)

  ! ################
  ! # Single terms #
  ! ################

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc =
  ! # yc =
  ! # zc =
  ! # vc = p
  ! # wc = rho
  ! # tc = T
  ! # sc = vis
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc =
  ! # pc =
  ! #################################

  CALL REDUCE( imax, jmax, kmax, u,   nstatavg, statavg, rc )
  CALL REDUCE( imax, jmax, kmax, v,   nstatavg, statavg, qc )
  CALL REDUCE( imax, jmax, kmax, w,   nstatavg, statavg, oc )
  CALL REDUCE( imax, jmax, kmax, p,   nstatavg, statavg, vc )
  CALL REDUCE( imax, jmax, kmax, rho, nstatavg, statavg, wc )
  CALL REDUCE( imax, jmax, kmax, T,   nstatavg, statavg, tc )

  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_U(j)   = MA_U(j)  + wrk2d(j,1)
    MA_V(j)   = MA_V(j)  + wrk2d(j,2)
    MA_W(j)   = MA_W(j)  + wrk2d(j,3)
    MA_P(j)   = MA_P(j)  + wrk2d(j,4)
    MA_R(j)   = MA_R(j)  + wrk2d(j,5)
    MA_T(j)   = MA_T(j)  + wrk2d(j,6)
  ENDDO

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
    CALL REDUCE( imax, jmax, kmax, vis, nstatavg, statavg, sc )
    CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,7), wrk2d(1,11) )
    DO j = 1,NNstat
      MA_VIS(j)  = MA_VIS(j) + wrk2d(j,7)
    ENDDO
  ELSE
    DO j = 1,NNstat
      MA_VIS(j)  = MA_VIS(j) + M_REAL(g(3)%size)
    ENDDO
  ENDIF

  ! ######################
  ! # Double cross-terms #
  ! ######################

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc  = u*u
  ! # yc  = v*v
  ! # zc  = w*w
  ! # vc  = p*p
  ! # wc  = rho*rho
  ! # tc  = T*T
  ! # sc  = vis*vis
  ! # rc  = u
  ! # qc  = v
  ! # oc = w
  ! # uc = rho*T*T
  ! #################################

  DO j = 1,NNstat*kmax
    uc(j,1,1) = wc(j,1,1)*tc(j,1,1)*tc(j,1,1)
    xc(j,1,1) = rc(j,1,1)*rc(j,1,1)
    yc(j,1,1) = qc(j,1,1)*qc(j,1,1)
    zc(j,1,1) = oc(j,1,1)*oc(j,1,1)
    vc(j,1,1) = vc(j,1,1)*vc(j,1,1)
    wc(j,1,1) = wc(j,1,1)*wc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,8), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UU(j) = MA_UU(j) + wrk2d(j,1)
    MA_VV(j) = MA_VV(j) + wrk2d(j,2)
    MA_WW(j) = MA_WW(j) + wrk2d(j,3)
    MA_PP(j) = MA_PP(j) + wrk2d(j,4)
    MA_RR(j) = MA_RR(j) + wrk2d(j,5)
    MA_TT(j) = MA_TT(j) + wrk2d(j,6)
    MA_RTT(j) = MA_RTT(j) + wrk2d(j,8)
  ENDDO

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
    DO j = 1,NNstat*kmax
      sc(j,1,1) = sc(j,1,1)*sc(j,1,1)
    ENDDO
    CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,7), wrk2d(1,11) )
    DO j = 1,NNstat
      MA_VIS2(j) = MA_VIS2(j) + wrk2d(j,7)
    ENDDO
  ELSE
    DO j = 1,NNstat
      MA_VIS2(j) = MA_VIS2(j) + M_REAL(g(3)%size)
    ENDDO
  ENDIF


  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc = u*v
  ! # yc = u*w
  ! # zc = v*w
  ! # vc = rho
  ! # wc = rho*rho
  ! # tc = rho*p
  ! # sc = rho*T
  ! # rc = rho*u
  ! # qc = rho*v
  ! # oc = rho*w
  ! # uc =
  ! # pc =
  ! #################################

  DO j = 1,NNstat*kmax
    xc(j,1,1) = rc(j,1,1)*qc(j,1,1)
    yc(j,1,1) = rc(j,1,1)*oc(j,1,1)
    zc(j,1,1) = qc(j,1,1)*oc(j,1,1)
  ENDDO

  CALL REDUCE( imax, jmax, kmax, rho, nstatavg, statavg, vc )
  CALL REDUCE( imax, jmax, kmax, p,   nstatavg, statavg, tc )
  CALL REDUCE( imax, jmax, kmax, T,   nstatavg, statavg, sc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = vc(j,1,1)*tc(j,1,1)
    sc(j,1,1) = vc(j,1,1)*sc(j,1,1)
    rc(j,1,1) = vc(j,1,1)*rc(j,1,1)
    qc(j,1,1) = vc(j,1,1)*qc(j,1,1)
    oc(j,1,1) = vc(j,1,1)*oc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,7), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,8), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UV(j) = MA_UV(j) + wrk2d(j,1)
    MA_UW(j) = MA_UW(j) + wrk2d(j,2)
    MA_VW(j) = MA_VW(j) + wrk2d(j,3)
    MA_RU(j) = MA_RU(j) + wrk2d(j,6)
    MA_RV(j) = MA_RV(j) + wrk2d(j,7)
    MA_RW(j) = MA_RW(j) + wrk2d(j,8)
    MA_RP(j) = MA_RP(j) + wrk2d(j,4)
    MA_RT(j) = MA_RT(j) + wrk2d(j,5)
  ENDDO

  ! ######################
  ! # Triple cross-terms #
  ! ######################

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc = rho*u*v
  ! # yc = rho*u*w
  ! # zc = rho*v*w
  ! # vc = rho
  ! # wc = u
  ! # tc = v
  ! # sc = w
  ! # rc = rho*u*u
  ! # qc = rho*v*v
  ! # oc = rho*w*w
  ! # uc =
  ! # pc =
  ! #################################

  DO j = 1,NNstat*kmax
    xc(j,1,1) = vc(j,1,1)*xc(j,1,1)
    yc(j,1,1) = vc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = vc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL REDUCE( imax, jmax, kmax, u, nstatavg, statavg, wc )
  CALL REDUCE( imax, jmax, kmax, v, nstatavg, statavg, tc )
  CALL REDUCE( imax, jmax, kmax, w, nstatavg, statavg, sc )

  DO j = 1,NNstat*kmax
    rc(j,1,1) = rc(j,1,1)*wc(j,1,1)
    qc(j,1,1) = qc(j,1,1)*tc(j,1,1)
    oc(j,1,1) = oc(j,1,1)*sc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUU(j) = MA_RUU(j) + wrk2d(j,1)
    MA_RVV(j) = MA_RVV(j) + wrk2d(j,2)
    MA_RWW(j) = MA_RWW(j) + wrk2d(j,3)
    MA_RUV(j) = MA_RUV(j) + wrk2d(j,4)
    MA_RUW(j) = MA_RUW(j) + wrk2d(j,5)
    MA_RVW(j) = MA_RVW(j) + wrk2d(j,6)
  ENDDO

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc = p
  ! # yc = T
  ! # zc = rho*v*w
  ! # vc = rho
  ! # wc = p*u
  ! # tc = p*v
  ! # sc = p*w
  ! # rc = T*u
  ! # qc = T*v
  ! # oc = T*w
  ! # uc =
  ! # pc =
  ! #################################

  CALL REDUCE( imax, jmax, kmax, p, nstatavg, statavg, xc )
  CALL REDUCE( imax, jmax, kmax, T, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    rc(j,1,1) = yc(j,1,1)*wc(j,1,1)
    qc(j,1,1) = yc(j,1,1)*tc(j,1,1)
    oc(j,1,1) = yc(j,1,1)*sc(j,1,1)
    wc(j,1,1) = xc(j,1,1)*wc(j,1,1)
    tc(j,1,1) = xc(j,1,1)*tc(j,1,1)
    sc(j,1,1) = xc(j,1,1)*sc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PU(j) = MA_PU(j) + wrk2d(j,1)
    MA_PV(j) = MA_PV(j) + wrk2d(j,2)
    MA_PW(j) = MA_PW(j) + wrk2d(j,3)
    MA_TU(j) = MA_TU(j) + wrk2d(j,4)
    MA_TV(j) = MA_TV(j) + wrk2d(j,5)
    MA_TW(j) = MA_TW(j) + wrk2d(j,6)
  ENDDO

  ! ############
  ! # Clean-up #
  ! ############

  ! ############################################################
  ! #                      DERIVATIVES U                       #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = du/dx*du/dx
  ! # wc = du/dy*du/dy
  ! # tc = du/dz*du/dz
  ! # sc = field du/dx
  ! # rc = field du/dy
  ! # qc = field du/dz
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, sc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, rc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, qc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, sc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, rc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, qc, nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*xc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*yc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Ux(j) = MA_Ux(j) + wrk2d(j,1)
    MA_Uy(j) = MA_Uy(j) + wrk2d(j,2)
    MA_Uz(j) = MA_Uz(j) + wrk2d(j,3)
    MA_Ux2(j) = MA_Ux2(j) + wrk2d(j,4)
    MA_Uy2(j) = MA_Uy2(j) + wrk2d(j,5)
    MA_Uz2(j) = MA_Uz2(j) + wrk2d(j,6)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = p*du/dx & u*p*du/dx
  ! # wc = p*du/dy & v*p*du/dx
  ! # tc = p*du/dz & p*(u*du/dx+v*du/dy+w*du/dz) & w*p*du/dx
  ! # sc = field du/dx
  ! # rc = field du/dy
  ! # qc = field du/dz
  ! # oc = u
  ! # uc = v
  ! # pc = w
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, p, nstatavg, statavg, tc )
  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, oc )
  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, uc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, pc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*tc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*tc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PUx(j) = MA_PUx(j) + wrk2d(j,1)
    MA_PUy(j) = MA_PUy(j) + wrk2d(j,2)
    MA_PUz(j) = MA_PUz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = C_2_R*oc(j,1,1)*vc(j,1,1)+uc(j,1,1)*wc(j,1,1)+pc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PHI1(j) = MA_PHI1(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = vc(j,1,1)*pc(j,1,1)
    wc(j,1,1) = vc(j,1,1)*uc(j,1,1)
    vc(j,1,1) = vc(j,1,1)*oc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PHI2(j) = MA_PHI2(j) + wrk2d(j,1)
    MA_PHI3(j) = MA_PHI3(j) + wrk2d(j,2)
    MA_PHI4(j) = MA_PHI4(j) + wrk2d(j,2)
    MA_PHI5(j) = MA_PHI5(j) + wrk2d(j,3)
    MA_PHI6(j) = MA_PHI6(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = (u & v & w)*du/dx
  ! # wc = (u & v & w)*du/dy
  ! # tc = (u & v & w)*du/dz
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc = u
  ! # uc = v
  ! # pc = w
  ! ######################################################

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*oc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*oc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*oc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UUX(j) = MA_UUX(j) + wrk2d(j,1)
    MA_UUY(j) = MA_UUY(j) + wrk2d(j,2)
    MA_UUZ(j) = MA_UUZ(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*uc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*uc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*uc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_VUX(j) = MA_VUX(j) + wrk2d(j,1)
    MA_VUY(j) = MA_VUY(j) + wrk2d(j,2)
    MA_VUZ(j) = MA_VUZ(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*pc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*pc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*pc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_WUX(j) = MA_WUX(j) + wrk2d(j,1)
    MA_WUY(j) = MA_WUY(j) + wrk2d(j,2)
    MA_WUZ(j) = MA_WUZ(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = rho*(u & v & w)*du/dx
  ! # wc = rho*(u & v & w)*du/dy
  ! # tc = rho*(u & v & w)*du/dz
  ! # sc = field du/dx
  ! # rc = field du/dy
  ! # qc = field du/dz
  ! # oc = rho*u
  ! # uc = rho*v
  ! # pc = rho*w
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = oc(j,1,1)*tc(j,1,1)
    uc(j,1,1) = uc(j,1,1)*tc(j,1,1)
    pc(j,1,1) = pc(j,1,1)*tc(j,1,1)
  ENDDO

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*oc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*oc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*oc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUUx(j) = MA_RUUx(j) + C_2_R*wrk2d(j,1)
    MA_RUUy(j) = MA_RUUy(j) + C_2_R*wrk2d(j,2)
    MA_RUUz(j) = MA_RUUz(j) + C_2_R*wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*uc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*uc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*uc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVx(j) = MA_RUVx(j) + wrk2d(j,1)
    MA_RUVy(j) = MA_RUVy(j) + wrk2d(j,2)
    MA_RUVz(j) = MA_RUVz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*pc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*pc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*pc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUWx(j) = MA_RUWx(j) + wrk2d(j,1)
    MA_RUWy(j) = MA_RUWy(j) + wrk2d(j,2)
    MA_RUWz(j) = MA_RUWz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = rho*u*(u & v & w)*du/dx & rho*v*v*du/dx
  ! # wc = rho*v*(u & v & w)*du/dy & rho*w*w*du/dx
  ! # tc = rho*w*(u & v & w)*du/dz & rho*v*w*du/dx
  ! # sc = field du/dx
  ! # rc = field du/dy
  ! # qc = field du/dz
  ! # oc = rho*u
  ! # uc = rho*v
  ! # pc = rho*w
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*oc(j,1,1)*tc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*uc(j,1,1)*tc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*pc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUUUkk(j) = MA_RUUUkk(j) +C_3_R*wrk2d(j,1) +C_2_R*wrk2d(j,2) +C_2_R*wrk2d(j,3)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*oc(j,1,1)*tc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*uc(j,1,1)*tc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*pc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVUkk(j) = MA_RUVUkk(j) + C_2_R*wrk2d(j,1) + wrk2d(j,2) + wrk2d(j,3)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*oc(j,1,1)*tc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*uc(j,1,1)*tc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*pc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUWUkk(j) = MA_RUWUkk(j) + C_2_R*wrk2d(j,1) + wrk2d(j,2) +wrk2d(j,3)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*uc(j,1,1)*tc(j,1,1)
  ENDDO
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    wc(j,1,1) = xc(j,1,1)*pc(j,1,1)*tc(j,1,1)
    tc(j,1,1) = xc(j,1,1)*uc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVVUkk(j) = MA_RVVUkk(j) + wrk2d(j,1)
    MA_RWWUkk(j) = MA_RWWUkk(j) + wrk2d(j,2)
    MA_RVWUkk(j) = MA_RVWUkk(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = rho*du/dx
  ! # wc = rho*du/dy
  ! # tc = rho*du/dz
  ! # sc = field du/dx
  ! # rc = field du/dy
  ! # qc = field du/dz
  ! # oc = rho*u
  ! # uc = rho*v
  ! # pc = rho*w
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*tc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*tc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUx(j) = MA_RUx(j) + wrk2d(j,1)
    MA_RUy(j) = MA_RUy(j) + wrk2d(j,2)
    MA_RUz(j) = MA_RUz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc =
  ! # wc =
  ! # tc =
  ! # sc = field tau_xx
  ! # rc = field tau_yy
  ! # qc = field tau_zz
  ! # oc = field tau_xy
  ! # uc = field tau_xz
  ! # pc = field tau_yz
  ! ######################################################

  DO j = 1,imax*jmax*kmax
    oc(j,1,1) =       rc(j,1,1)
    uc(j,1,1) =       qc(j,1,1)
    rc(j,1,1) =      -sc(j,1,1)
    qc(j,1,1) =      -sc(j,1,1)
    sc(j,1,1) = C_2_R*sc(j,1,1)
  ENDDO

  ! ############################################################
  ! #                      DERIVATIVES V                       #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = du/dy*dv/dx & du/dx*dv/dy
  ! # sc = field tau_xx
  ! # rc = field tau_yy
  ! # qc = field tau_zz
  ! # oc = field tau_xy
  ! # uc = field tau_xz
  ! # pc = field tau_yz
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, wc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tc, wrk3d, wrk2d,wrk3d)

  DO j = 1,imax*jmax*kmax
    oc(j,1,1) = oc(j,1,1) +       wc(j,1,1)
    sc(j,1,1) = sc(j,1,1) -       tc(j,1,1)
    rc(j,1,1) = rc(j,1,1) + C_2_R*tc(j,1,1)
    qc(j,1,1) = qc(j,1,1) -       tc(j,1,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, wc, nstatavg, statavg, vc )
  CALL REDUCE( imax,jmax,kmax, tc, nstatavg, statavg, wc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc,  wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc,  wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Vx(j) = MA_Vx(j) + wrk2d(j,1)
    MA_Vy(j) = MA_Vy(j) + wrk2d(j,2)
    MA_VxUy(j) = MA_VxUy(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = xc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_UXVY(j) = MA_UXVY(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = p
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = p*dv/dx & u*p*dv/dx & p*dv/dy & (u & v & w)*p*dv/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, p, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_PVX(j) = MA_PVX(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = tc(j,1,1)*yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_PHI4(j) = MA_PHI4(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_PVY(j) = MA_PVY(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = tc(j,1,1)*yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_PHI2(j) = MA_PHI2(j) + wrk2d(j,1)
    MA_PHI1(j) = MA_PHI1(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = tc(j,1,1)*yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_PHI3(j) = MA_PHI3(j) +       wrk2d(j,1)
    MA_PHI4(j) = MA_PHI4(j) + C_2_R*wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = tc(j,1,1)*yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_PHI5(j) = MA_PHI5(j) + wrk2d(j,1)
    MA_PHI6(j) = MA_PHI6(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = rho
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = rho*dv/dx & rho*dv/dy & z1*rho*dv/dx & z1*rho*dv/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVx(j) = MA_RVx(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVy(j) = MA_RVy(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = u & rho*u
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = (u & rho*u)*(dv/dx & dv/dy)
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UVX(j) = MA_UVX(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_UVY(j) = MA_UVY(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    yc(j,1,1) = yc(j,1,1)*tc(j,1,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVx(j) = MA_RUVx(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUVy(j) = MA_RUVy(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = rho*u
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = rho*u*(u & v & w)*(dv/dx & dv/dy)
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUVUkk(j) = MA_RUVUkk(j) + wrk2d(j,1)
  ENDDO
  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUUUkk(j) = MA_RUUUkk(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVVUkk(j) = MA_RVVUkk(j) + C_2_R*wrk2d(j,1)
  ENDDO
  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUVUkk(j) = MA_RUVUkk(j) + C_2_R*wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWUkk(j) = MA_RVWUkk(j) + wrk2d(j,1)
  ENDDO
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUWUkk(j) = MA_RUWUkk(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = v & rho*v
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = (v & rho*v)*(dv/dx & dv/dy )
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_VVx(j) = MA_VVx(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_VVy(j) = MA_VVy(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    yc(j,1,1) = yc(j,1,1)*tc(j,1,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVVx(j) = MA_RVVx(j) + C_2_R*wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVVy(j) = MA_RVVy(j) + C_2_R*wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = rho*v
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = rho*v*(v & w )*dv/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVVUkk(j) = MA_RVVUkk(j) + C_3_R*wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWUkk(j) = MA_RVWUkk(j) + C_2_R*wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = w & rho*w
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = (w & rho*w)*(dv/dx & dv/dy )
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_WVx(j) = MA_WVx(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_WVy(j) = MA_WVy(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    yc(j,1,1) = yc(j,1,1)*tc(j,1,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWx(j) = MA_RVWx(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWy(j) = MA_RVWy(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = rho*w
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = rho*w*(w )*dv/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RWWUkk(j) = MA_RWWUkk(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = dv/dx*dv/dx
  ! # zc = du/dz
  ! # vc = dv/dx
  ! # wc = dv/dy
  ! # tc = dv/dy*dv/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  DO j = 1,NNstat*kmax
    yc(j,1,1) = vc(j,1,1)*vc(j,1,1)
    tc(j,1,1) = wc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, yc,  wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,2), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_Vx2(j) = MA_Vx2(j) + wrk2d(j,1)
    MA_Vy2(j) = MA_Vy2(j) + wrk2d(j,2)
  ENDDO

  ! ############################################################
  ! #                      DERIVATIVES W                       #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = dw/dz
  ! # zc = du/dz
  ! # vc = dw/dx
  ! # wc = dv/dy
  ! # tc = dw/dx*du/dz & dw/dz*du/dx & dw/dz*dv/dy
  ! # sc = field tau_xx
  ! # rc = field tau_yy
  ! # qc = field tau_zz
  ! # oc = field tau_xy
  ! # uc = field tau_xz
  ! # pc = field tau_yz
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, vc, wrk3d, wrk2d,wrk3d)

  DO j = 1,imax*jmax*kmax
    uc(j,1,1) = uc(j,1,1) +       tc(j,1,1)
    sc(j,1,1) = sc(j,1,1) -       vc(j,1,1)
    rc(j,1,1) = rc(j,1,1) -       vc(j,1,1)
    qc(j,1,1) = qc(j,1,1) + C_2_R*vc(j,1,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, tc, nstatavg, statavg, vc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = vc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Wz(j) = MA_Wz(j) + wrk2d(j,1)
    MA_Wx(j) = MA_Wx(j) + wrk2d(j,2)
    MA_WxUz(j) = MA_WxUz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = yc(j,1,1)*xc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_UxWz(j) = MA_UxWz(j) + wrk2d(j,1)
  ENDDO

  DO j = 1,NNstat*kmax
    tc(j,1,1) = wc(j,1,1)*yc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_VyWz(j) = MA_VyWz(j) + wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = u*p*dw/dz
  ! # yc = dw/dz
  ! # zc = p*dw/dz
  ! # vc = dw/dx
  ! # wc = p*dw/dx & u*p*dw/dx & v*p*dw/dz
  ! # tc = w*p*dw/dz
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, p,   nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = zc(j,1,1)*vc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*yc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PWz(j) = MA_PWz(j) + wrk2d(j,1)
    MA_PWx(j) = MA_PWx(j) + wrk2d(j,2)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, xc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = wc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PHI6(j) = MA_PHI6(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, wc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    xc(j,1,1) = xc(j,1,1)*zc(j,1,1)
    wc(j,1,1) = wc(j,1,1)*zc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PHI2(j) = MA_PHI2(j) +       wrk2d(j,1)
    MA_PHI1(j) = MA_PHI1(j) +       wrk2d(j,1)
    MA_PHI3(j) = MA_PHI3(j) +       wrk2d(j,2)
    MA_PHI4(j) = MA_PHI4(j) +       wrk2d(j,2)
    MA_PHI5(j) = MA_PHI5(j) +       wrk2d(j,3)
    MA_PHI6(j) = MA_PHI6(j) + C_2_R*wrk2d(j,3)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = rho*dw/dz
  ! # yc = dw/dz
  ! # zc =
  ! # vc = dw/dx
  ! # wc =
  ! # tc = rho*dw/dx
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, xc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = xc(j,1,1)*vc(j,1,1)
    xc(j,1,1) = xc(j,1,1)*yc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RWx(j) = MA_RWx(j) + wrk2d(j,1)
    MA_RWz(j) = MA_RWz(j) + wrk2d(j,2)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = rho
  ! # yc = dw/dz
  ! # zc = u*dw/dz & rho*u*dw/dz
  ! # vc = dw/dx
  ! # wc = u*dw/dx & rho*u*dw/dx
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, xc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = zc(j,1,1)*vc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*yc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UWz(j) = MA_UWz(j) + wrk2d(j,1)
    MA_UWx(j) = MA_UWx(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    wc(j,1,1) = xc(j,1,1)*wc(j,1,1)
    zc(j,1,1) = xc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUWz(j) = MA_RUWz(j) + wrk2d(j,1)
    MA_RUWx(j) = MA_RUWx(j) + wrk2d(j,2)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = (u & v & w )*rho*u*dw/dz
  ! # yc = dw/dz
  ! # zc = rho*u*dw/dz
  ! # vc = dw/dx
  ! # wc = rho*u*dw/dx
  ! # tc = (u & v & w)*rho*u*dw/dx
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = tc(j,1,1)*zc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUUUkk(j) = MA_RUUUkk(j) + wrk2d(j,1)
    MA_RUWUkk(j) = MA_RUWUkk(j) + wrk2d(j,2)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = tc(j,1,1)*zc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUVUkk(j) = MA_RUVUkk(j) + wrk2d(j,1)
    MA_RVWUkk(j) = MA_RVWUkk(j) + wrk2d(j,2)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = tc(j,1,1)*zc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*wc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RUWUkk(j) = MA_RUWUkk(j) + C_2_R*wrk2d(j,1)
    MA_RWWUkk(j) = MA_RWWUkk(j) + C_2_R*wrk2d(j,2)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = rho
  ! # yc = dw/dz
  ! # zc = v*dw/dz & rho*v*dw/dz
  ! # vc = dw/dx
  ! # wc = v*dw/dx & rho*v*dw/dx
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, xc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = zc(j,1,1)*vc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*yc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_VWz(j) = MA_VWz(j) + wrk2d(j,1)
    MA_VWx(j) = MA_VWx(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    wc(j,1,1) = xc(j,1,1)*wc(j,1,1)
    zc(j,1,1) = xc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVWz(j) = MA_RVWz(j) + wrk2d(j,1)
    MA_RVWx(j) = MA_RVWx(j) + wrk2d(j,2)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = (v & w )*rho*v*dw/dz
  ! # yc = dw/dz
  ! # zc = rho*v*dw/dz
  ! # vc = dw/dx
  ! # wc = rho*v*dw/dx
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = tc(j,1,1)*zc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVVUkk(j) = MA_RVVUkk(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = tc(j,1,1)*zc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWUkk(j) = MA_RVWUkk(j) + C_2_R*wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = rho
  ! # yc = dw/dz
  ! # zc = w*dw/dz & rho*w*dw/dz
  ! # vc = dw/dx
  ! # wc = w*dw/dx & rho*w*dw/dx
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, w,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, xc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = zc(j,1,1)*vc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*yc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_WWz(j) = MA_WWz(j) + wrk2d(j,1)
    MA_WWx(j) = MA_WWx(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    wc(j,1,1) = wc(j,1,1)*xc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RWWz(j) = MA_RWWz(j) + C_2_R*wrk2d(j,1)
    MA_RWWx(j) = MA_RWWx(j) + C_2_R*wrk2d(j,2)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = rho*w*(w )*dw/dz
  ! # yc = dw/dz
  ! # zc = rho*w*dw/dz
  ! # vc = dw/dx
  ! # wc = rho*w*dw/dx
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = tc(j,1,1)*zc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RWWUkk(j) = MA_RWWUkk(j) + C_3_R*wrk2d(j,1)
  ENDDO

  ! ###################################################
  ! # Temporary array storage
  ! #
  ! # xc = dw/dz*dw/dz
  ! # yc = dw/dz
  ! # zc = dw/dx*dw/dx
  ! # vc = dw/dx
  ! # wc =
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  DO j = 1,NNstat*kmax
    xc(j,1,1) = yc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = vc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, xc,  wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc,  wrk2d(1,2), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_Wz2(j) = MA_Wz2(j) + wrk2d(j,1)
    MA_Wx2(j) = MA_Wx2(j) + wrk2d(j,2)
  ENDDO

  ! ############################################################
  ! #                   CLOSING DERIVATIVES V-W                #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = dv/dz*dw/dy
  ! # vc = dv/dz*dv/dz
  ! # wc = dw/dy*dw/dy
  ! # tc =
  ! # sc = field tau_xx
  ! # rc = field tau_yy
  ! # qc = field tau_zz
  ! # oc = field tau_xy
  ! # uc = field tau_xz
  ! # pc = field tau_yz
  ! ######################################################

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, zc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, yc, wrk3d, wrk2d,wrk3d)

  DO j = 1,imax*jmax*kmax
    pc(j,1,1) = zc(j,1,1) + yc(j,1,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, yc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, zc, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    zc(j,1,1) = xc(j,1,1)*yc(j,1,1)
    vc(j,1,1) = xc(j,1,1)*xc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*yc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,5), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Vz(j) = MA_Vz(j) + wrk2d(j,1)
    MA_Wy(j) = MA_Wy(j) + wrk2d(j,2)
    MA_WyVz(j) = MA_WyVz(j) + wrk2d(j,3)
    MA_Vz2(j) = MA_Vz2(j) + wrk2d(j,4)
    MA_Wy2(j) = MA_Wy2(j) + wrk2d(j,5)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc =
  ! # vc =
  ! # wc = rho*dv/dz
  ! # tc = rho*dw/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, wc )

  DO j = 1,NNstat*kmax
    tc(j,1,1) = wc(j,1,1)*yc(j,1,1)
    wc(j,1,1) = wc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVz(j) = MA_RVz(j) + wrk2d(j,1)
    MA_RWy(j) = MA_RWy(j) + wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = p*dv/dz & w*p*dv/dz
  ! # vc = p*dw/dy & v*p*dw/dy
  ! # wc = v
  ! # tc = w
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, p, nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = zc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PVZ(j) = MA_PVZ(j) + wrk2d(j,1)
    MA_PWY(j) = MA_PWY(j) + wrk2d(j,2)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, wc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = vc(j,1,1)*wc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PHI4(j) = MA_PHI4(j) + wrk2d(j,1)
    MA_PHI6(j) = MA_PHI6(j) + wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = u*dv/dz & rho*u*dv/dz
  ! # vc = u*dw/dy & rho*u*dw/dy
  ! # wc = rho
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, wc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = zc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UVZ(j) = MA_UVZ(j) + wrk2d(j,1)
    MA_UWY(j) = MA_UWY(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    zc(j,1,1) = zc(j,1,1)*wc(j,1,1)
    vc(j,1,1) = vc(j,1,1)*wc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVz(j) = MA_RUVz(j) + wrk2d(j,1)
    MA_RUWy(j) = MA_RUWy(j) + wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = rho*u*dv/dz
  ! # vc = rho*u*dw/dy
  ! # wc = w*rho*u*dv/dz
  ! # tc = v*rho*u*dw/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, wc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = wc(j,1,1)*zc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*vc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVUkk(j) = MA_RUVUkk(j) + wrk2d(j,1)
    MA_RUWUkk(j) = MA_RUWUkk(j) + wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = v*dv/dz & rho*v*dv/dz
  ! # vc = v*dw/dy & rho*v*dw/dy
  ! # wc = rho
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, wc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = zc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_VVZ(j) = MA_VVZ(j) + wrk2d(j,1)
    MA_VWY(j) = MA_VWY(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    zc(j,1,1) = zc(j,1,1)*wc(j,1,1)
    vc(j,1,1) = vc(j,1,1)*wc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVVz(j)  = MA_RVVz(j)  + C_2_R*wrk2d(j,1)
    MA_RVWy(j) = MA_RVWy(j) +       wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = rho*v*dv/dz
  ! # vc = rho*v*dw/dy
  ! # wc = w*rho*v*dv/dz
  ! # tc = (v & w )*rho*v*dw/dy
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    tc(j,1,1) = tc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWUkk(j) = MA_RVWUkk(j) + wrk2d(j,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    wc(j,1,1) = tc(j,1,1)*zc(j,1,1)
    tc(j,1,1) = tc(j,1,1)*vc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,2), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVVUkk(j) = MA_RVVUkk(j) + C_2_R*wrk2d(j,1)
    MA_RWWUkk(j) = MA_RWWUkk(j) + C_2_R*wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = w*dv/dz & rho*w*dv/dz
  ! # vc = w*dw/dy & rho*w*dw/dy
  ! # wc = rho
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, w,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, wc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = zc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*xc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_WVZ(j) = MA_WVZ(j) + wrk2d(j,1)
    MA_WWY(j) = MA_WWY(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    zc(j,1,1) = zc(j,1,1)*wc(j,1,1)
    vc(j,1,1) = vc(j,1,1)*wc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVWz(j) = MA_RVWz(j) +       wrk2d(j,1)
    MA_RWWy(j)  = MA_RWWy(j)  + C_2_R*wrk2d(j,2)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dz
  ! # yc = dw/dy
  ! # zc = rho*w*dv/dz
  ! # vc = rho*w*dw/dy
  ! # wc = (w )*rho*w*dv/dz
  ! # tc =
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )
  DO j = 1,NNstat*kmax
    wc(j,1,1) = tc(j,1,1)*zc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,1), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_RVWUkk(j) = MA_RVWUkk(j) + wrk2d(j,1)
  ENDDO

  ! ############################################################
  ! #              STRESS TENSOR AND ITS DERIVATIVES           #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc =
  ! # yc =
  ! # zc =
  ! # vc =
  ! # wc =
  ! # tc =
  ! # sc = field tau_xx
  ! # rc = field tau_yy
  ! # qc = field tau_zz
  ! # oc = field tau_xy
  ! # uc = field tau_xz
  ! # pc = field tau_yz
  ! ######################################################

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
    DO j = 1,imax*jmax*kmax
      sc(j,1,1) = vis(j,1,1)*sc(j,1,1)*c23
      rc(j,1,1) = vis(j,1,1)*rc(j,1,1)*c23
      qc(j,1,1) = vis(j,1,1)*qc(j,1,1)*c23
      oc(j,1,1) = vis(j,1,1)*oc(j,1,1)
      uc(j,1,1) = vis(j,1,1)*uc(j,1,1)
      pc(j,1,1) = vis(j,1,1)*pc(j,1,1)
    ENDDO
  ELSE
    DO j = 1,imax*jmax*kmax
      sc(j,1,1) = sc(j,1,1)*c23
      rc(j,1,1) = rc(j,1,1)*c23
      qc(j,1,1) = qc(j,1,1)*c23
    ENDDO
  ENDIF

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = d(tau_xk)/dk
  ! # yc = d(tau_yk)/dk
  ! # zc = d(tau_zk)/dk
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), sc, vc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), oc, wc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), uc, tc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, wc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, tc, nstatavg, statavg, zc )

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), oc, vc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, wc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = xc(j,1,1) + wc(j,1,1)
  ENDDO

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), rc, vc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, wc )
  DO j = 1,NNstat*kmax
    yc(j,1,1) = yc(j,1,1) + wc(j,1,1)
  ENDDO

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), pc, vc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, wc )
  DO j = 1,NNstat*kmax
    zc(j,1,1) = zc(j,1,1) + wc(j,1,1)
  ENDDO

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), uc, vc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, wc )
  DO j = 1,NNstat*kmax
    xc(j,1,1) = xc(j,1,1) + wc(j,1,1)
  ENDDO

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), pc, vc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, wc )
  DO j = 1,NNstat*kmax
    yc(j,1,1) = yc(j,1,1) + wc(j,1,1)
  ENDDO

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), qc, vc, wrk3d, wrk2d,wrk3d)
  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, wc )
  DO j = 1,NNstat*kmax
    zc(j,1,1) = zc(j,1,1) + wc(j,1,1)
  ENDDO

  DO j = 1,NNstat*kmax
    xc(j,1,1) = visc*xc(j,1,1)
    yc(j,1,1) = visc*yc(j,1,1)
    zc(j,1,1) = visc*zc(j,1,1)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, sc, nstatavg, statavg, vc )
  CALL REDUCE( imax,jmax,kmax, rc, nstatavg, statavg, wc )
  CALL REDUCE( imax,jmax,kmax, qc, nstatavg, statavg, tc )
  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, sc )
  CALL REDUCE( imax,jmax,kmax, uc, nstatavg, statavg, rc )
  CALL REDUCE( imax,jmax,kmax, pc, nstatavg, statavg, qc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = visc*vc(j,1,1)
    wc(j,1,1) = visc*wc(j,1,1)
    tc(j,1,1) = visc*tc(j,1,1)
    sc(j,1,1) = visc*sc(j,1,1)
    rc(j,1,1) = visc*rc(j,1,1)
    qc(j,1,1) = visc*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,7), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,8), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,9), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUXkk(j) = MA_TAUXkk(j)  + wrk2d(j,1)
    MA_TAUYkk(j) = MA_TAUYkk(j)  + wrk2d(j,2)
    MA_TAUZkk(j) = MA_TAUZkk(j)  + wrk2d(j,3)
    MA_TAUxx(j) = MA_TAUxx(j) + wrk2d(j,4)
    MA_TAUyy(j) = MA_TAUyy(j) + wrk2d(j,5)
    MA_TAUzz(j) = MA_TAUzz(j) + wrk2d(j,6)
    MA_TAUxy(j) = MA_TAUxy(j) + wrk2d(j,7)
    MA_TAUxz(j) = MA_TAUxz(j) + wrk2d(j,8)
    MA_TAUyz(j) = MA_TAUyz(j) + wrk2d(j,9)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = d(tau_xk)/dk
  ! # yc = d(tau_yk)/dk
  ! # zc = d(tau_zk)/dk
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc = (u & v & w)*d(tau_xk)/dk
  ! # uc = (u & v & w)*d(tau_yk)/dk
  ! # pc = (u & v & w)*d(tau_zk)/dk
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, pc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = pc(j,1,1)*xc(j,1,1)
    uc(j,1,1) = pc(j,1,1)*yc(j,1,1)
    pc(j,1,1) = pc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UTAUXkk(j) = MA_UTAUXkk(j) + wrk2d(j,1)
    MA_UTAUYkk(j) = MA_UTAUYkk(j) + wrk2d(j,2)
    MA_UTAUZkk(j) = MA_UTAUZkk(j) + wrk2d(j,3)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, pc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = pc(j,1,1)*xc(j,1,1)
    uc(j,1,1) = pc(j,1,1)*yc(j,1,1)
    pc(j,1,1) = pc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_VTAUXkk(j) = MA_VTAUXkk(j) + wrk2d(j,1)
    MA_VTAUYkk(j) = MA_VTAUYkk(j) + wrk2d(j,2)
    MA_VTAUZkk(j) = MA_VTAUZkk(j) + wrk2d(j,3)
  ENDDO

  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, pc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = pc(j,1,1)*xc(j,1,1)
    uc(j,1,1) = pc(j,1,1)*yc(j,1,1)
    pc(j,1,1) = pc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_WTAUXkk(j) = MA_WTAUXkk(j) + wrk2d(j,1)
    MA_WTAUYkk(j) = MA_WTAUYkk(j) + wrk2d(j,2)
    MA_WTAUZkk(j) = MA_WTAUZkk(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = d(tau_xk)/dk
  ! # yc = d(tau_yk)/dk
  ! # zc = d(tau_zk)/dk
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc = T*d(tau_xk)/dk
  ! # uc = T*d(tau_yk)/dk
  ! # pc = T*d(tau_zk)/dk
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, T, nstatavg, statavg, pc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = pc(j,1,1)*xc(j,1,1)
    uc(j,1,1) = pc(j,1,1)*yc(j,1,1)
    pc(j,1,1) = pc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TTAUXkk(j) = MA_TTAUXkk(j) + wrk2d(j,1)
    MA_TTAUYkk(j) = MA_TTAUYkk(j) + wrk2d(j,2)
    MA_TTAUZkk(j) = MA_TTAUZkk(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = du/dx
  ! # yc = du/dy
  ! # zc = du/dz
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc = tau_xk*du/dk
  ! # uc = tau_yk*du/dk & tau_zk*du/dk
  ! # pc = phi
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, oc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, uc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, pc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, uc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, pc, nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*xc(j,1,1) + sc(j,1,1)*yc(j,1,1) +           rc(j,1,1)*zc(j,1,1)
    uc(j,1,1) = sc(j,1,1)*xc(j,1,1) + wc(j,1,1)*yc(j,1,1) +           qc(j,1,1)*zc(j,1,1)
    pc(j,1,1) = oc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUXkUk(j) = MA_TAUXkUk(j) + wrk2d(j,1)
    MA_TAUYkUk(j) = MA_TAUYkUk(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    uc(j,1,1) = rc(j,1,1)*xc(j,1,1) + qc(j,1,1)*yc(j,1,1) +           tc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUZkUk(j) = MA_TAUZkUk(j) + wrk2d(j,1)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dv/dx
  ! # yc = dv/dy
  ! # zc = dv/dz
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc = tau_xk*dv/dk
  ! # uc = tau_yk*dv/dk & tau_zk*dv/dk
  ! # pc = phi
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, oc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, uc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, uc, nstatavg, statavg, yc )

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, oc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*xc(j,1,1) + sc(j,1,1)*yc(j,1,1) +           rc(j,1,1)*zc(j,1,1)
    uc(j,1,1) = sc(j,1,1)*xc(j,1,1) + wc(j,1,1)*yc(j,1,1) +           qc(j,1,1)*zc(j,1,1)
    pc(j,1,1) = pc(j,1,1) + uc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUXkVk(j) = MA_TAUXkVk(j) + wrk2d(j,1)
    MA_TAUYkVk(j) = MA_TAUYkVk(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    uc(j,1,1) = rc(j,1,1)*xc(j,1,1) + qc(j,1,1)*yc(j,1,1) +           tc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUZkVk(j) = MA_TAUZkVk(j) + wrk2d(j,1)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dw/dx
  ! # yc = dw/dy
  ! # zc = dw/dz
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc = tau_xk*dw/dk
  ! # uc = tau_yk*dw/dk & tau_zk*dw/dk
  ! # pc = phi
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, oc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, uc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, uc, nstatavg, statavg, yc )

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, oc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*xc(j,1,1) + sc(j,1,1)*yc(j,1,1) +           rc(j,1,1)*zc(j,1,1)
    uc(j,1,1) = sc(j,1,1)*xc(j,1,1) + wc(j,1,1)*yc(j,1,1) +           qc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUXkWk(j) = MA_TAUXkWk(j) + wrk2d(j,1)
    MA_TAUYkWk(j) = MA_TAUYkWk(j) + wrk2d(j,2)
  ENDDO

  DO j = 1,NNstat*kmax
    uc(j,1,1) = rc(j,1,1)*xc(j,1,1) + qc(j,1,1)*yc(j,1,1) +           tc(j,1,1)*zc(j,1,1)
    pc(j,1,1) = pc(j,1,1) + uc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,1), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TAUZkWk(j) = MA_TAUZkWk(j) + wrk2d(j,1)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = u*phi
  ! # yc = v*phi
  ! # zc = w*phi
  ! # vc = tau_xx
  ! # wc = tau_yy
  ! # tc = tau_zz
  ! # sc = tau_xy
  ! # rc = tau_xz
  ! # qc = tau_yz
  ! # oc =
  ! # uc =
  ! # pc = phi
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, zc )

  DO j = 1,NNstat*kmax
    xc(j,1,1) = pc(j,1,1)*xc(j,1,1)
    yc(j,1,1) = pc(j,1,1)*yc(j,1,1)
    zc(j,1,1) = pc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,4), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UPHI(j) = MA_UPHI(j) + wrk2d(j,1)
    MA_VPHI(j) = MA_VPHI(j) + wrk2d(j,2)
    MA_WPHI(j) = MA_WPHI(j) + wrk2d(j,3)
  ENDDO

  ! ############
  ! # Clean-up #
  ! ############

  ! ############################################################
  ! #                        DERIVATIVES RHO                   #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = untouched
  ! # wc = untouched
  ! # tc = untouched
  ! # sc = field drho/dx
  ! # rc = field drho/dy
  ! # qc = field drho/dz
  ! # oc =
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), rho, sc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), rho, rc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), rho, qc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, sc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, rc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, qc, nstatavg, statavg, zc )

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Rx(j) = MA_Rx(j) + wrk2d(j,1)
    MA_Ry(j) = MA_Ry(j) + wrk2d(j,2)
    MA_Rz(j) = MA_Rz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = u
  ! # wc = v
  ! # tc = w
  ! # sc = u*drho/dx
  ! # rc = u*drho/dy
  ! # qc = u*drho/dz
  ! # oc = (u & v & w )*u*drho/dx
  ! # uc = (u & v & w )*u*drho/dy
  ! # pc = (u & v & w )*u*drho/dz
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, vc )
  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, wc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    sc(j,1,1) = vc(j,1,1)*xc(j,1,1)
    rc(j,1,1) = vc(j,1,1)*yc(j,1,1)
    qc(j,1,1) = vc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_URx(j) = MA_URx(j) + wrk2d(j,1)
    MA_URy(j) = MA_URy(j) + wrk2d(j,2)
    MA_URz(j) = MA_URz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = vc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = vc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUUx(j) = MA_RUUx(j) + wrk2d(j,1)
    MA_RUUy(j) = MA_RUUy(j) + wrk2d(j,2)
    MA_RUUz(j) = MA_RUUz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = wc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = wc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVx(j) = MA_RUVx(j) + wrk2d(j,1)
    MA_RUVy(j) = MA_RUVy(j) + wrk2d(j,2)
    MA_RUVz(j) = MA_RUVz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = tc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = tc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUWx(j) = MA_RUWx(j) + wrk2d(j,1)
    MA_RUWy(j) = MA_RUWy(j) + wrk2d(j,2)
    MA_RUWz(j) = MA_RUWz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = u
  ! # wc = v
  ! # tc = w
  ! # sc = u*drho/dx
  ! # rc = u*drho/dy
  ! # qc = u*drho/dz
  ! # oc = u*u*u*drho/dx & v*u*u*drho/dx
  ! # uc = v*u*u*drho/dy & v*v*u*drho/dy
  ! # pc = w*u*u*drho/dz & v*w*u*drho/dz
  ! ######################################################

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*vc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*vc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*vc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUUUkk(j) = MA_RUUUkk(j) + wrk2d(j,1) +           wrk2d(j,2) +           wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*wc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*wc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*wc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUVUkk(j) = MA_RUVUkk(j) + wrk2d(j,1) +           wrk2d(j,2) +           wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = u
  ! # wc = v
  ! # tc = w
  ! # sc = v*drho/dx
  ! # rc = v*drho/dy
  ! # qc = v*drho/dz
  ! # oc = (v & w )*v*drho/dx
  ! # uc = (v & w )*v*drho/dy
  ! # pc = (v & w )*v*drho/dz
  ! ######################################################

  DO j = 1,NNstat*kmax
    sc(j,1,1) = wc(j,1,1)*xc(j,1,1)
    rc(j,1,1) = wc(j,1,1)*yc(j,1,1)
    qc(j,1,1) = wc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_VRx(j) = MA_VRx(j) + wrk2d(j,1)
    MA_VRy(j) = MA_VRy(j) + wrk2d(j,2)
    MA_VRz(j) = MA_VRz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = wc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = wc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVVx(j) = MA_RVVx(j) + wrk2d(j,1)
    MA_RVVy(j) = MA_RVVy(j) + wrk2d(j,2)
    MA_RVVz(j) = MA_RVVz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = tc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = tc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVWx(j)  = MA_RVWx(j)  + wrk2d(j,1)
    MA_RVWy(j) = MA_RVWy(j) + wrk2d(j,2)
    MA_RVWz(j) = MA_RVWz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = u
  ! # wc = v
  ! # tc = w
  ! # sc = v*drho/dx
  ! # rc = v*drho/dy
  ! # qc = v*drho/dz
  ! # oc = u*v*v*drho/dx & u*w*v*drho/dx
  ! # uc = v*v*v*drho/dy & v*w*v*drho/dy
  ! # pc = w*v*v*drho/dz & w*w*v*drho/dz
  ! ######################################################

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*wc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*wc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*wc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVVUkk(j) = MA_RVVUkk(j) + wrk2d(j,1) +           wrk2d(j,2) +           wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*tc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*tc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*tc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RVWUkk(j) = MA_RVWUkk(j) + wrk2d(j,1) +           wrk2d(j,2) +           wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = u
  ! # wc = v
  ! # tc = w
  ! # sc = w*drho/dx
  ! # rc = w*drho/dy
  ! # qc = w*drho/dz
  ! # oc = (w)*w*drho/dx
  ! # uc = (w)*w*drho/dy
  ! # pc = (w)*w*drho/dz
  ! ######################################################

  DO j = 1,NNstat*kmax
    sc(j,1,1) = tc(j,1,1)*xc(j,1,1)
    rc(j,1,1) = tc(j,1,1)*yc(j,1,1)
    qc(j,1,1) = tc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_WRx(j) = MA_WRx(j) + wrk2d(j,1)
    MA_WRy(j) = MA_WRy(j) + wrk2d(j,2)
    MA_WRz(j) = MA_WRz(j) + wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = tc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = tc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RWWx(j) = MA_RWWx(j) + wrk2d(j,1)
    MA_RWWy(j) = MA_RWWy(j) + wrk2d(j,2)
    MA_RWWz(j) = MA_RWWz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = drho/dx
  ! # yc = drho/dy
  ! # zc = drho/dz
  ! # vc = u
  ! # wc = v
  ! # tc = w
  ! # sc = w*drho/dx
  ! # rc = w*drho/dy
  ! # qc = w*drho/dz
  ! # oc = u*w*w*drho/dx & u*u*w*drho/dx
  ! # uc = v*w*w*drho/dy & u*v*w*drho/dy
  ! # pc = w*w*w*drho/dz & u*w*w*drho/dz
  ! ######################################################

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*tc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*tc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*tc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RWWUkk(j) = MA_RWWUkk(j) + wrk2d(j,1) +           wrk2d(j,2) +           wrk2d(j,3)
  ENDDO

  DO j = 1,NNstat*kmax
    oc(j,1,1) = vc(j,1,1)*vc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = wc(j,1,1)*vc(j,1,1)*rc(j,1,1)
    pc(j,1,1) = tc(j,1,1)*vc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUWUkk(j) = MA_RUWUkk(j) + wrk2d(j,1) +           wrk2d(j,2) +           wrk2d(j,3)
  ENDDO

  ! ############
  ! # Clean-up #
  ! ############

  ! ############################################################
  ! #                        DERIVATIVES P                     #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dp/dx
  ! # yc = dp/dy
  ! # zc = dp/dz
  ! # vc = u*dp/dx
  ! # wc = v*dp/dy
  ! # tc = w*dp/dz
  ! # sc = u*dp/dx + v*dp/dy + w*dp/dz
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), p, sc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), p, rc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), p, qc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, sc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, rc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, qc, nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, u, nstatavg, statavg, rc )
  CALL REDUCE( imax,jmax,kmax, v, nstatavg, statavg, qc )
  CALL REDUCE( imax,jmax,kmax, w, nstatavg, statavg, oc )

  DO j = 1,NNstat*kmax
    vc(j,1,1) = xc(j,1,1)*rc(j,1,1)
    wc(j,1,1) = yc(j,1,1)*qc(j,1,1)
    tc(j,1,1) = zc(j,1,1)*oc(j,1,1)
    sc(j,1,1) = vc(j,1,1) + wc(j,1,1) + tc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc,  wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc,  wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc,  wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc,  wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc,  wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,6), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, sc,  wrk2d(1,7), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Px(j)  = MA_Px(j)  + wrk2d(j,1)
    MA_Py(j)  = MA_Py(j)  + wrk2d(j,2)
    MA_Pz(j)  = MA_Pz(j)  + wrk2d(j,3)
    MA_UPx(j) = MA_UPx(j) + wrk2d(j,4)
    MA_VPy(j) = MA_VPy(j) + wrk2d(j,5)
    MA_WPz(j) = MA_WPz(j) + wrk2d(j,6)
    MA_UkPk(j) = MA_UkPk(j)  + wrk2d(j,7)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dp/dx
  ! # yc = dp/dy
  ! # zc = dp/dz
  ! # vc = u*dp/dx
  ! # wc = v*dp/dy
  ! # tc = w*dp/dz
  ! # sc = u*(u*dp/dx + v*dp/dy + w*dp/dz)
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc = v*(u*dp/dx + v*dp/dy + w*dp/dz)
  ! # pc = w*(u*dp/dx + v*dp/dy + w*dp/dz)
  ! ######################################################

  DO j = 1,NNstat*kmax
    pc(j,1,1) = oc(j,1,1)*sc(j,1,1)
    uc(j,1,1) = qc(j,1,1)*sc(j,1,1)
    sc(j,1,1) = rc(j,1,1)*sc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PHI1(j)  = MA_PHI1(j) + wrk2d(j,1)
    MA_PHI4(j)  = MA_PHI4(j) + wrk2d(j,2)
    MA_PHI6(j)  = MA_PHI6(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dp/dx
  ! # yc = dp/dy
  ! # zc = dp/dz
  ! # vc = u*dp/dy
  ! # wc = u*dp/dz
  ! # tc = v*dp/dx
  ! # sc = v*dp/dz
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc = w*dp/dx
  ! # pc = w*dp/dy
  ! ######################################################

  DO j = 1,NNstat*kmax
    vc(j,1,1) = yc(j,1,1)*rc(j,1,1)
    wc(j,1,1) = zc(j,1,1)*rc(j,1,1)
    tc(j,1,1) = xc(j,1,1)*qc(j,1,1)
    sc(j,1,1) = zc(j,1,1)*qc(j,1,1)
    uc(j,1,1) = xc(j,1,1)*oc(j,1,1)
    pc(j,1,1) = yc(j,1,1)*oc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_UPy(j) = MA_UPy(j) + wrk2d(j,1)
    MA_UPz(j) = MA_UPz(j) + wrk2d(j,2)
    MA_VPx(j) = MA_VPx(j) + wrk2d(j,3)
    MA_VPz(j) = MA_VPz(j) + wrk2d(j,4)
    MA_WPx(j) = MA_WPx(j) + wrk2d(j,5)
    MA_WPy(j) = MA_WPy(j) + wrk2d(j,6)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dp/dx
  ! # yc = dp/dy
  ! # zc = dp/dz
  ! # vc = T
  ! # wc = T*dp/dx
  ! # tc = T*dp/dy
  ! # sc = T*dp/dz
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, T, nstatavg, statavg, vc )

  DO j = 1,NNstat*kmax
    wc(j,1,1) = vc(j,1,1)*xc(j,1,1)
    tc(j,1,1) = vc(j,1,1)*yc(j,1,1)
    sc(j,1,1) = vc(j,1,1)*zc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, wc,  wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc,  wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, sc,  wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_TPx(j) = MA_TPx(j) + wrk2d(j,1)
    MA_TPy(j) = MA_TPy(j) + wrk2d(j,2)
    MA_TPz(j) = MA_TPz(j) + wrk2d(j,3)
  ENDDO

  ! ############################################################
  ! #                        DERIVATIVES T                     #
  ! ############################################################

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = dT/dx
  ! # yc = dT/dy
  ! # zc = dT/dz
  ! # vc = field dT/dx
  ! # wc = field dT/dy
  ! # tc = field dT/dz
  ! # sc = untouched
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), T, vc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), T, wc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), T, tc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, wc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, tc, nstatavg, statavg, zc )

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Tx(j)  = MA_Tx(j)  + wrk2d(j,1)
    MA_Ty(j)  = MA_Ty(j)  + wrk2d(j,2)
    MA_Tz(j)  = MA_Tz(j)  + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = p*dT/dx
  ! # yc = p*dT/dy
  ! # zc = p*dT/dz
  ! # vc = field dT/dx
  ! # wc = field dT/dy
  ! # tc = field dT/dz
  ! # sc = p
  ! # rc = u
  ! # qc = v
  ! # oc = w
  ! # uc =
  ! ######################################################

  CALL REDUCE( imax,jmax,kmax, p, nstatavg, statavg, sc )

  DO j = 1,NNstat*kmax
    xc(j,1,1) = xc(j,1,1)*sc(j,1,1)
    yc(j,1,1) = yc(j,1,1)*sc(j,1,1)
    zc(j,1,1) = zc(j,1,1)*sc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_PTx(j) = MA_PTx(j) + wrk2d(j,1)
    MA_PTy(j) = MA_PTy(j) + wrk2d(j,2)
    MA_PTz(j) = MA_PTz(j) + wrk2d(j,3)
  ENDDO

  ! ####################################################
  ! # Temporary array storage
  ! #
  ! # xc = field d2T/dx2
  ! # yc = field d2T/dy2
  ! # zc = field d2T/dz2
  ! # vc = d2T/dx2
  ! # wc = d2T/dy2
  ! # tc = d2T/dz2
  ! # sc = d2T/dx2 + d2T/dy2 + d2T/dz2
  ! # rc = u*(d2T/dx2 + d2T/dy2 + d2T/dz2)
  ! # qc = v*(d2T/dx2 + d2T/dy2 + d2T/dz2)
  ! # oc = w*(d2T/dx2 + d2T/dy2 + d2T/dz2)
  ! # uc =
  ! # pc =
  ! ######################################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), vc, xc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), wc, yc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tc, zc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, xc, nstatavg, statavg, vc )
  CALL REDUCE( imax,jmax,kmax, yc, nstatavg, statavg, wc )
  CALL REDUCE( imax,jmax,kmax, zc, nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    sc(j,1,1)  =  vc(j,1,1) + wc(j,1,1) + tc(j,1,1)
    rc(j,1,1)  =  rc(j,1,1)*sc(j,1,1)
    qc(j,1,1)  =  qc(j,1,1)*sc(j,1,1)
    oc(j,1,1) = oc(j,1,1)*sc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, sc,  wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc,  wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc,  wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,4), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_Tkk(j)  = MA_Tkk(j)  + wrk2d(j,1)
    MA_UTkk(j) = MA_UTkk(j) + wrk2d(j,2)
    MA_VTkk(j) = MA_VTkk(j) + wrk2d(j,3)
    MA_WTkk(j) = MA_WTkk(j) + wrk2d(j,4)
  ENDDO

  ! ############
  ! # Clean-up #
  ! ############

  ! ############################################################
  ! #                     SKEWNESS & FLATNESS                  #
  ! ############################################################

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc = rho
  ! # yc = u
  ! # zc = v
  ! # vc = w
  ! # wc = p
  ! # tc = T
  ! # sc = rho^3 and rho^4
  ! # rc = u^3 and u^4
  ! # qc = v^3 and v^4
  ! # oc = w^3 and w^4
  ! # uc = p^3 and p^4
  ! # pc = T^3 and T^4
  ! #################################

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, u,   nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, v,   nstatavg, statavg, zc )
  CALL REDUCE( imax,jmax,kmax, w,   nstatavg, statavg, vc )
  CALL REDUCE( imax,jmax,kmax, p,   nstatavg, statavg, wc )
  CALL REDUCE( imax,jmax,kmax, T,   nstatavg, statavg, tc )

  DO j = 1,NNstat*kmax
    sc(j,1,1) = xc(j,1,1)**C_3_R
    rc(j,1,1) = yc(j,1,1)**C_3_R
    qc(j,1,1) = zc(j,1,1)**C_3_R
    oc(j,1,1) = vc(j,1,1)**C_3_R
    uc(j,1,1) = wc(j,1,1)**C_3_R
    pc(j,1,1) = tc(j,1,1)**C_3_R
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,6), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_R3(j) = MA_R3(j) + wrk2d(j,1)
    MA_U3(j) = MA_U3(j) + wrk2d(j,2)
    MA_V3(j) = MA_V3(j) + wrk2d(j,3)
    MA_W3(j) = MA_W3(j) + wrk2d(j,4)
    MA_P3(j) = MA_P3(j) + wrk2d(j,5)
    MA_T3(j) = MA_T3(j) + wrk2d(j,6)
  ENDDO

  DO j = 1,NNstat*kmax
    sc(j,1,1) = sc(j,1,1)*xc(j,1,1)
    rc(j,1,1) = rc(j,1,1)*yc(j,1,1)
    qc(j,1,1) = qc(j,1,1)*zc(j,1,1)
    oc(j,1,1) = oc(j,1,1)*vc(j,1,1)
    uc(j,1,1) = uc(j,1,1)*wc(j,1,1)
    pc(j,1,1) = pc(j,1,1)*tc(j,1,1)
  ENDDO
  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, oc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, uc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, pc, wrk2d(1,6), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_R4(j) = MA_R4(j) + wrk2d(j,1)
    MA_U4(j) = MA_U4(j) + wrk2d(j,2)
    MA_V4(j) = MA_V4(j) + wrk2d(j,3)
    MA_W4(j) = MA_W4(j) + wrk2d(j,4)
    MA_P4(j) = MA_P4(j) + wrk2d(j,5)
    MA_T4(j) = MA_T4(j) + wrk2d(j,6)
  ENDDO

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile, 'LEAVING DNS_SAVE_AVGIJ' )
#endif

  ! ############
  ! # Clean-up #
  ! ############

  ! ############################################################
  ! #         Additional terms for the T'2 equation            #
  ! ############################################################

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc = p
  ! # yc = T
  ! # zc =
  ! # vc = p*T*dudx
  ! # wc = p*T*dvdy
  ! # tc = p*T*dwdz
  ! # sc = T*dudx
  ! # rc = T*dvdy
  ! # qc = T*dwdz
  ! # oc = dudx
  ! # uc = dvdy
  ! # pc = dwdz
  ! #################################

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, qc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, oc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, uc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, uc, nstatavg, statavg, pc )
  CALL REDUCE( imax,jmax,kmax, oc, nstatavg, statavg, uc )
  CALL REDUCE( imax,jmax,kmax, qc, nstatavg, statavg, oc )
  CALL REDUCE( imax,jmax,kmax, p, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, T, nstatavg, statavg, yc )

  DO j = 1,NNstat*kmax
    sc(j,1,1) = yc(j,1,1)*oc(j,1,1)
    rc(j,1,1) = yc(j,1,1)*uc(j,1,1)
    qc(j,1,1) = yc(j,1,1)*pc(j,1,1)
    vc(j,1,1) = xc(j,1,1)*sc(j,1,1)
    wc(j,1,1) = xc(j,1,1)*rc(j,1,1)
    tc(j,1,1) = xc(j,1,1)*qc(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, sc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, rc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, qc, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, vc, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, wc, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tc, wrk2d(1,6), wrk2d(1,11) )
  DO j = 1,NNstat
    MA_TUx(j) = MA_TUx(j) + wrk2d(j,1)
    MA_TVy(j) = MA_TVy(j) + wrk2d(j,2)
    MA_TWz(j) = MA_TWz(j) + wrk2d(j,3)
    MA_PTUx(j) = MA_PTUx(j) + wrk2d(j,4)
    MA_PTVy(j) = MA_PTVy(j) + wrk2d(j,5)
    MA_PTWz(j) = MA_PTWz(j) + wrk2d(j,6)
  ENDDO

  ! #################################
  ! # Temporary array storage
  ! #
  ! # xc =
  ! # yc =
  ! # zc =
  ! # vc = d(rho*u*T*T)/dx
  ! # wc = d(rho*v*T*T)/dy
  ! # tc = d(rho*w*T*T)/dz
  ! # sc =
  ! # rc =
  ! # qc =
  ! # oc =
  ! # uc =
  ! # pc =
  ! #################################

  DO j = 1,imax*jmax*kmax
    xc(j,1,1) = rho(j,1,1)*u(j,1,1)*T(j,1,1)**C_2_R
    yc(j,1,1) = rho(j,1,1)*v(j,1,1)*T(j,1,1)**C_2_R
    zc(j,1,1) = rho(j,1,1)*w(j,1,1)*T(j,1,1)**C_2_R
  ENDDO
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), xc, vc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), yc, wc, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), zc, tc, wrk3d, wrk2d,wrk3d)

  CALL REDUCE( imax,jmax,kmax, vc, nstatavg, statavg, xc )
  CALL REDUCE( imax,jmax,kmax, wc, nstatavg, statavg, yc )
  CALL REDUCE( imax,jmax,kmax, tc, nstatavg, statavg, zc )

  CALL SUM1V1D_V( NNstat, kmax, xc, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, yc, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, zc, wrk2d(1,3), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_RUTTx(j) = MA_RUTTx(j) + wrk2d(j,1)
    MA_RVTTy(j) = MA_RVTTy(j) + wrk2d(j,2)
    MA_RWTTz(j) = MA_RWTTz(j) + wrk2d(j,3)
  ENDDO

  RETURN
END SUBROUTINE AVG_FLOW_ZT_REDUCE

! ###################################################################
! ###################################################################
SUBROUTINE AVG_TKE_ZT_REDUCE(rho, u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, mean1d, wrk2d)

  ! ##############################################
  ! # Running statistics for kinetic energy before
  ! # filtering (Favre averages)
  ! #
  ! # 10/12/2000 Juan Pedro Mellado
  ! ##############################################

  USE TLAB_VARS
  USE AVGS, ONLY: SUM1V1D_V

  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax) :: rho, u, v, w
  TREAL, DIMENSION(imax,jmax,kmax) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7

  TREAL mean1d(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)

  TINTEGER j
  TINTEGER NNstat

  NNstat = jmax*nstatavg

  ! #################################
  ! # Temporary array storage
  ! #
  ! # tmp1 = rho
  ! # tmp2 = rho*u*u
  ! # tmp3 = rho*v*v
  ! # tmp4 = rho*w*w
  ! # tmp5 = rho*u
  ! # tmp6 = rho*v
  ! # tmp7 = rho*w
  ! #################################

  CALL REDUCE( imax, jmax, kmax, rho, nstatavg, statavg, tmp1 )
  CALL REDUCE( imax, jmax, kmax, u,   nstatavg, statavg, tmp2 )
  CALL REDUCE( imax, jmax, kmax, v,   nstatavg, statavg, tmp3 )
  CALL REDUCE( imax, jmax, kmax, w,   nstatavg, statavg, tmp4 )

  DO j = 1,NNstat*kmax
    tmp5(j,1,1) = tmp1(j,1,1)*tmp2(j,1,1)
    tmp6(j,1,1) = tmp1(j,1,1)*tmp3(j,1,1)
    tmp7(j,1,1) = tmp1(j,1,1)*tmp4(j,1,1)
    tmp2(j,1,1) = tmp5(j,1,1)*tmp2(j,1,1)
    tmp3(j,1,1) = tmp6(j,1,1)*tmp3(j,1,1)
    tmp4(j,1,1) = tmp7(j,1,1)*tmp4(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, tmp5, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp6, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp7, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp2, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp3, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp4, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
    MA_FLT_RU(j) = MA_FLT_RU(j) + wrk2d(j,1)
    MA_FLT_RV(j) = MA_FLT_RV(j) + wrk2d(j,2)
    MA_FLT_RW(j) = MA_FLT_RW(j) + wrk2d(j,3)
    MA_FLT_RUU(j) = MA_FLT_RUU(j) + wrk2d(j,4)
    MA_FLT_RVV(j) = MA_FLT_RVV(j) + wrk2d(j,5)
    MA_FLT_RWW(j) = MA_FLT_RWW(j) + wrk2d(j,6)
  ENDDO

  RETURN
END SUBROUTINE AVG_TKE_ZT_REDUCE
