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
! #    Derivatives in z1
! #
! # 09/25/00 Juan Pedro Mellado
! #####################################################

MODULE AVG_SCAL_ZT

  USE TLAB_CONSTANTS, ONLY : efile
#ifdef TRACE_ON
  USE TLAB_CONSTANTS, ONLY : tfile
  USE TLAB_PROCS,     ONLY : TLAB_WRITE_ASCII
#endif
  USE TLAB_VARS, ONLY : isize_field, imax,jmax,kmax,inb_scal, isize_wrk2d, imode_eqns
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : nstatavg, statavg
  USE TLAB_VARS, ONLY : itransport, visc
  USE AVGS, ONLY: SUM1V1D_V

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: AVG_SCAL_ZT_REDUCE

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE(q,z1, hq,txc, mean1d_sc, wrk2d,wrk3d )
    IMPLICIT NONE

    TREAL, DIMENSION(isize_field,*), INTENT(IN   ), TARGET :: q, z1
    TREAL, DIMENSION(isize_field,*), INTENT(INOUT), TARGET :: txc, hq
    TREAL mean1d_sc(nstatavg,jmax,MS_SCALAR_SIZE,*)
    TREAL wrk2d(isize_wrk2d,*)
    TREAL wrk3d(*)

    ! -------------------------------------------------------------------
    TINTEGER is, bcs(2,1)
    TINTEGER NNstat

    ! Pointers to existing allocated space
    TREAL, DIMENSION(:), POINTER :: u, v, w, rho, p, vis
    TREAL, DIMENSION(:), POINTER :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12

    ! ###################################################################
#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'ENTERING AVG_SCAL_ZT_REDUCE' )
#endif

    bcs = 0
    NNstat = nstatavg*jmax

    ! Define pointers
    u   => q(:,1)
    v   => q(:,2)
    w   => q(:,3)
    IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
       rho => q(:,5)
       p   => q(:,6)
       IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)
    ENDIF

    tmp1 => hq(:,1)
    tmp2 => hq(:,2)
    tmp3 => hq(:,3)
    tmp4 => txc(:,1)
    tmp5 => txc(:,2)
    tmp6 => txc(:,3)
    tmp7 => txc(:,4)
    tmp8 => txc(:,5)
    tmp9 => txc(:,6)
    tmp10 => txc(:,7)
    tmp11 => txc(:,8)
    tmp12 => txc(:,9)

    ! ################################################################
    ! # Mean Terms
    ! ################################################################

#define m_rho tmp1
#define m_u   tmp2
#define m_v   tmp3
#define m_w   tmp4
#define m_z1  tmp6

    CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, m_rho)
    CALL REDUCE( imax,jmax,kmax, u,   nstatavg, statavg, m_u)
    CALL REDUCE( imax,jmax,kmax, v,   nstatavg, statavg, m_v)
    CALL REDUCE( imax,jmax,kmax, w,   nstatavg, statavg, m_w)

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_M1(NNstat, m_rho, m_u, m_v, m_w, &
          m_z1, z1(1,is), tmp7, tmp8, tmp9, tmp10, tmp11, &
          mean1d_sc(1,1,1,is), wrk2d)
    ENDDO

    ! ##################################################################
    ! # Density Gradient terms
    ! ##################################################################

#define m_rho_x       tmp9
#define m_rho_y       tmp10

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), rho, tmp7, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), rho, tmp8, wrk3d, wrk2d,wrk3d)

    CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_rho_x)
    CALL REDUCE(imax,jmax,kmax, tmp8, nstatavg, statavg, m_rho_y)

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_R1(NNstat, m_rho_x, m_rho_y, m_u, &
          m_v, m_w, m_z1, z1(1,is), tmp5, tmp7, tmp8, &
          tmp11, tmp12, wrk3d, mean1d_sc(1,1,1,is), wrk2d)
    ENDDO

#undef m_rho_x
#undef m_rho_y

    ! ################################################################
    ! # U - V gradient terms
    ! ################################################################

#define m_u_x        tmp9
#define m_v_y        tmp10

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp7, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp8, wrk3d, wrk2d,wrk3d)

    CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_u_x)
    CALL REDUCE(imax,jmax,kmax, tmp8, nstatavg, statavg, m_v_y)

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_UV1(NNstat, m_u_x, m_v_y, m_rho, &
          m_u, m_v, m_w, m_z1, z1(1,is), tmp5, tmp7, &
          tmp8, tmp11, tmp12, wrk3d, mean1d_sc(1,1,1,is),   wrk2d)
    ENDDO

#undef m_u_x
#undef m_v_y

    ! ################################################################
    ! # W gradients terms
    ! ################################################################

#define m_w_x        tmp9
#define m_w_y        tmp10

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tmp7, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp8, wrk3d, wrk2d,wrk3d)

    CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_w_x)
    CALL REDUCE(imax,jmax,kmax, tmp8, nstatavg, statavg, m_w_y)

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_W1(NNstat, m_w_x, m_w_y, m_rho, &
          m_u, m_v, m_z1, z1(1,is), tmp7, &
          tmp8, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk2d)
    ENDDO

#undef m_w_x
#undef m_w_y

    ! ###############################################################
    ! # U -V second gradient parts
    ! ###############################################################

#define m_v_x        tmp9
#define m_u_y        tmp10

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp7, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp8, wrk3d, wrk2d,wrk3d)

    CALL REDUCE(imax,jmax,kmax, tmp7, nstatavg, statavg, m_v_x)
    CALL REDUCE(imax,jmax,kmax, tmp8, nstatavg, statavg, m_u_y)

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_UV2(NNstat, m_v_x, m_u_y, m_rho, &
          m_u, m_v, m_z1, z1(1,is), tmp7, &
          tmp8, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk2d)
    ENDDO

#undef m_v_x
#undef m_u_y

    ! ####################################################################
    ! # Scalar Gradients and Dissipation
    ! ####################################################################

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_S2(NNstat, m_rho, m_u, m_v, m_w, m_z1, p, vis, z1(1,is), &
          tmp5, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk2d, wrk3d)
    ENDDO

    ! ####################################################################
    ! # Viscous stresses terms
    ! ####################################################################

#undef m_rho
#undef m_u
#undef m_v
#undef m_w

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_TS1(NNstat, m_z1, u, v, w, vis, z1(1,is), tmp1, tmp2, tmp3, &
          tmp4, tmp5, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk2d, wrk3d)
    ENDDO

    ! ##################################################################
    ! #  Pressure gradient terms
    ! ##################################################################

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), p, tmp7, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), p, tmp8, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), p, tmp9, wrk3d, wrk2d,wrk3d)

#define m_p_x tmp10
#define m_p_y tmp11
#define m_p_z tmp12

    CALL REDUCE( imax,jmax,kmax, tmp7,  nstatavg, statavg, m_p_x)
    CALL REDUCE( imax,jmax,kmax, tmp8,  nstatavg, statavg, m_p_y)
    CALL REDUCE( imax,jmax,kmax, tmp9,  nstatavg, statavg, m_p_z)

    DO is=1, inb_scal
      CALL AVG_SCAL_ZT_REDUCE_P1(NNstat, m_p_x, m_p_y, m_p_z, &
          m_z1, z1(1,is), tmp7, tmp8, tmp9, &
          mean1d_sc(1,1,1,is), wrk2d)
    ENDDO

#undef m_p_x
#undef m_p_y
#undef m_p_z

#ifdef TRACE_ON
    CALL TLAB_WRITE_ASCII(tfile, 'LEAVING AVG_SCAL_ZT_REDUCE' )
#endif

    RETURN
  END SUBROUTINE AVG_SCAL_ZT_REDUCE

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_M1(NNstat, m_rho, m_u, m_v, m_w, &
      m_z1, z1, tmp1, tmp2, tmp3, tmp4, tmp5, mean1d_sc, wrk2d)
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

    CALL SUM1V1D_V( NNstat, kmax, m_z1, wrk2d(1,1), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1, wrk2d(1,2), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u, wrk2d(1,3), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v, wrk2d(1,4), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_w, wrk2d(1,5), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_z1, wrk2d(1,6), wrk2d(1,11))

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_2, wrk2d(1,1), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v_2, wrk2d(1,2), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_w_2, wrk2d(1,3), wrk2d(1,11))

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_v, wrk2d(1,1), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_w, wrk2d(1,2), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v_w, wrk2d(1,3), wrk2d(1,11))

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_z1_2, wrk2d(1,1), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_v_z1_2, wrk2d(1,2), wrk2d(1,11))
    CALL SUM1V1D_V( NNstat, kmax, m_rho_w_z1_2, wrk2d(1,3), wrk2d(1,11))

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

    CALL SUM1V1D_V( NNstat, kmax, m_z1_z1, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_z1_u, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_z1_v, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_z1_w, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_z1_gamma, wrk2d(1,1), wrk2d(1,11) )

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
    CALL SUM1V1D_V( NNstat, kmax, m_z1_3, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_z1_4, wrk2d(1,2), wrk2d(1,11) )
    DO j = 1,NNstat
      MS_S3(j) = MS_S3(j) + wrk2d(j,1)
      MS_S4(j) = MS_S4(j) + wrk2d(j,2)
    ENDDO

#undef m_z1_3
#undef m_z1_4

#undef m_z1

    RETURN
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_M1

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_P1(NNstat, m_p_x, m_p_y, m_p_z, &
      m_z1, z1, m_p_x_z1, m_p_y_z1, m_p_z_z1, mean1d_sc, wrk2d)
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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_P1

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_R1(NNstat, m_rho_x, m_rho_y, &
      m_u, m_v, m_w, m_z1, z1, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, mean1d_sc, wrk2d)
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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_z1, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_z1, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_u, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_v, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_2_u, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_2_v, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_v, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_u, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_z1_w, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_z1_w, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_u_u_z1, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_v_v_z1, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_u_v_z1, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_u_v_z1, wrk2d(1,4), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_x_u_w_z1, wrk2d(1,5), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_y_v_w_z1, wrk2d(1,6), wrk2d(1,11) )

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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_R1

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_S2(NNstat, m_rho, m_u, m_v, m_w, m_z1, p, vis, z1, &
      tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, mean1d_sc, wrk2d, wrk3d)
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
    TREAL tmp6(*)
    TREAL tmp7(*)

    TREAL p(*)
    TREAL vis(*)
    TREAL z1(imax,jmax,kmax)
    TREAL mean1d_sc(nstatavg,jmax,*)
    TREAL wrk2d(isize_wrk2d,*)
    TREAL wrk3d(*)

    TINTEGER j, bcs(2,1)

    bcs = 0

    CALL REDUCE( imax,jmax,kmax, z1, nstatavg, statavg, m_z1)

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1, tmp2, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1, tmp3, wrk3d, wrk2d,wrk3d)
    CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), z1, tmp4, wrk3d, wrk2d,wrk3d)

    IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
      DO j=1, kmax*imax*jmax
        tmp5(j) = vis(j)*tmp2(j)
      ENDDO
    ELSE
      DO j=1, kmax*imax*jmax
        tmp5(j) = tmp2(j)
      ENDDO
    ENDIF

    CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp5, tmp6, wrk3d, wrk2d,wrk3d)

#define m_fxx tmp7

    CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_fxx)

    CALL SUM1V1D_V( NNstat, kmax, m_fxx, wrk2d(1,1),         wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_z1_fxx, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_u_fxx, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_v_fxx, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_w_fxx, wrk2d(1,4), wrk2d(1,11) )

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
    CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5, tmp6, wrk3d, wrk2d,wrk3d)

#define m_fyy tmp7

    CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_fyy)

    CALL SUM1V1D_V( NNstat, kmax, m_fyy, wrk2d(1,1),         wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_z1_fyy, wrk2d(1,1),         wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_u_fyy, wrk2d(1,2),         wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_v_fyy, wrk2d(1,3),         wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_w_fyy, wrk2d(1,4),         wrk2d(1,11) )

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

    CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp5, tmp6, wrk3d, wrk2d,wrk3d)

#define m_fzz tmp7

    CALL REDUCE(imax,jmax,kmax, tmp6, nstatavg, statavg, m_fzz)

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

    CALL SUM1V1D_V( NNstat, kmax, m_z1_fzz, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_u_fzz,  wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_v_fzz,  wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_w_fzz,  wrk2d(1,4), wrk2d(1,11) )

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

    CALL REDUCE(imax,jmax,kmax, tmp2, nstatavg, statavg, m_z1_x)
    CALL REDUCE(imax,jmax,kmax, tmp3, nstatavg, statavg, m_z1_y)
    CALL REDUCE(imax,jmax,kmax, tmp4, nstatavg, statavg, m_z1_z)

#define m_vis  tmp1
    IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
      CALL REDUCE(imax,jmax,kmax, vis, nstatavg, statavg, m_vis)
    ELSE
      DO j = 1,NNstat*kmax
        m_vis(j) = C_1_R
      ENDDO
    ENDIF

    CALL SUM1V1D_V( NNstat, kmax, m_z1_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_z1_y, wrk2d(1,2), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_vis_z1_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_vis_z1_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_vis_z1_z, wrk2d(1,3), wrk2d(1,11) )

    DO j = 1,NNstat
      MS_Fx(j) = MS_Fx(j) + wrk2d(j,1)
      MS_Fy(j) = MS_Fy(j) + wrk2d(j,2)
      MS_Fz(j) = MS_Fz(j) + wrk2d(j,3)
    ENDDO

    ! here m_vis_z1_z contains the total dissipation
#define m_eps m_vis_z1_z
    DO j = 1,NNstat*kmax
      m_eps(j) = m_vis(j)*(m_z1_x(j)**2  +m_z1_y(j)**2+m_z1_z(j)**2)
      m_vis_z1_x(j) = m_vis_z1_x(j)*m_z1(j)
      m_vis_z1_y(j) = m_vis_z1_y(j)*m_z1(j)
    ENDDO

    CALL SUM1V1D_V( NNstat, kmax, m_vis_z1_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_vis_z1_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_eps,      wrk2d(1,3), wrk2d(1,11) )

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
    CALL REDUCE( imax,jmax,kmax, p,  nstatavg, statavg, m_p)

#define m_p_z1_x tmp3
#define m_p_z1_y wrk3d
#define m_p_z1_z tmp2

    DO j = 1,NNstat*kmax
      m_p_z1_x(j) = m_p(j)*m_z1_x(j)
      m_p_z1_y(j) = m_p(j)*m_z1_y(j)
      m_p_z1_z(j) = m_p(j)*m_z1_z(j)
    ENDDO

    CALL SUM1V1D_V( NNstat, kmax, m_p_z1_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_p_z1_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_p_z1_z, wrk2d(1,3), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_x,    wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_y,    wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_z1_x, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_z1_y, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_x_u,    wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_y_v,    wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_z1_x_u, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_z1_y_v, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_x_v, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_y_u, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_x_w, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_y_w, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_u_z1_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_v_v_z1_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_v_z1_x, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_v_z1_y, wrk2d(1,4), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_w_z1_x, wrk2d(1,5), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_v_w_z1_y, wrk2d(1,6), wrk2d(1,11) )

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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_S2

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_TS1(NNstat, m_z1, u, v, w, vis, z1, tmp1, tmp2, tmp3, tmp4, tmp5, &
      tmp6, tmp7, tmp8, tmp9, tmp10, tmp11, mean1d_sc, wrk2d, wrk3d)
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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_TS1

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_UV1(NNstat, m_u_x, m_v_y, &
      m_rho, m_u, m_v, m_w, m_z1, z1, tmp1, tmp2, tmp3, tmp4, &
      tmp5, tmp6, mean1d_sc, wrk2d)
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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_2_u_x, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_2_v_y, wrk2d(1,4), wrk2d(1,11) )

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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_u_x_z1, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_v_v_y_z1, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_x_v_z1, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_v_y_z1, wrk2d(1,4), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_x_w_z1, wrk2d(1,5), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_v_y_w_z1, wrk2d(1,6), wrk2d(1,11) )

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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_UV1

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_UV2(NNstat, m_v_x, m_u_y, &
      m_rho, m_u, m_v, m_z1, z1, tmp1, tmp2, tmp3, tmp4, &
      mean1d_sc, wrk2d)
    IMPLICIT NONE

    TINTEGER NNstat
    TREAL m_v_x(*)
    TREAL m_u_y(*)
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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_v_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_u_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_v_x_z1, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_y_v_z1, wrk2d(1,4), wrk2d(1,11) )

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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_UV2

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE AVG_SCAL_ZT_REDUCE_W1(NNstat, m_w_x, m_w_y, &
      m_rho, m_u, m_v, m_z1, z1, tmp1, tmp2, tmp3, &
      tmp4, mean1d_sc, wrk2d)
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

    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_w_x, wrk2d(1,1), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_z1_w_y, wrk2d(1,2), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_u_w_x_z1, wrk2d(1,3), wrk2d(1,11) )
    CALL SUM1V1D_V( NNstat, kmax, m_rho_v_w_y_z1, wrk2d(1,4), wrk2d(1,11) )

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
  END SUBROUTINE AVG_SCAL_ZT_REDUCE_W1

END MODULE AVG_SCAL_ZT
