#include "types.h"
#include "dns_const.h"
#include "avgij_map.h"

SUBROUTINE DNS_SAVE_SCBDGIJ(rho,u,v,w,p,z1, vis, &
     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12, &
     mean1d_sc, wrk1d,wrk2d,wrk3d )

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

  USE DNS_GLOBAL

  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax)   :: u, v, w, p, rho, vis
  TREAL, DIMENSION(imax,jmax,kmax,*) :: z1
  TREAL, DIMENSION(*)                :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(*)                :: tmp7,tmp8,tmp9,tmp10,tmp11,tmp12

  TREAL mean1d_sc(nstatavg,jmax,MS_SCALAR_SIZE,*)

  TREAL wrk1d(*)
  TREAL wrk2d(isize_wrk2d,*)
  TREAL wrk3d(*)

  TINTEGER is, bcs(2,1)
  TINTEGER NNstat

  bcs = 0
  NNstat = nstatavg*jmax

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING DNS_SAVE_SCBDGIJ' )
#endif

! ################################################################
! # Mean Terms
! ################################################################

#define m_rho tmp1
#define m_u   tmp2
#define m_v   tmp3
#define m_w   tmp4
#define m_z1  tmp6

  CALL REDUCE( imax,jmax,kmax, rho, nstatavg, statavg, m_rho)
  CALL REDUCE( imax,jmax,kmax, u,  nstatavg, statavg, m_u)
  CALL REDUCE( imax,jmax,kmax, v,  nstatavg, statavg, m_v)
  CALL REDUCE( imax,jmax,kmax, w,  nstatavg, statavg, m_w)

  DO is=1, inb_scal
     CALL DNS_SAVE_SCBDGIJ_M1(NNstat, m_rho, m_u, m_v, m_w, &
          m_z1, z1(1,1,1,is), tmp7, tmp8, tmp9, tmp10, tmp11, &
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
     CALL DNS_SAVE_SCBDGIJ_R1(NNstat, m_rho_x, m_rho_y, m_u, &
          m_v, m_w, m_z1, z1(1,1,1,is), tmp5, tmp7, tmp8, &
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
     CALL DNS_SAVE_SCBDGIJ_UV1(NNstat, m_u_x, m_v_y, m_rho, &
          m_u, m_v, m_w, m_z1, z1(1,1,1,is), tmp5, tmp7, &
          tmp8, tmp11, tmp12, wrk3d, mean1d_sc(1,1,1,is), &
          wrk2d)
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
     CALL DNS_SAVE_SCBDGIJ_W1(NNstat, m_w_x, m_w_y, m_rho, &
          m_u, m_v, m_w, m_z1, z1(1,1,1,is), tmp7, &
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
     CALL DNS_SAVE_SCBDGIJ_UV2(NNstat, m_v_x, m_u_y, m_rho, &
          m_u, m_v, m_w, m_z1, z1(1,1,1,is), tmp7, &
          tmp8, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk2d)
  ENDDO

#undef m_v_x
#undef m_u_y

! ####################################################################
! # Scalar Gradients and Dissipation
! ####################################################################

  DO is=1, inb_scal
     CALL DNS_SAVE_SCBDGIJ_S2(NNstat, m_rho, m_u, m_v, m_w, m_z1, p, vis, z1(1,1,1,is), &
          tmp5, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk1d, wrk2d, wrk3d)
  ENDDO

! ####################################################################
! # Viscous stresses terms
! ####################################################################

#undef m_rho
#undef m_u
#undef m_v
#undef m_w

  DO is=1, inb_scal
     CALL DNS_SAVE_SCBDGIJ_TS1(NNstat, m_z1, u, v, w, vis, z1(1,1,1,is), tmp1, tmp2, tmp3, &
          tmp4, tmp5, tmp7, tmp8, tmp9, tmp10, tmp11, tmp12, mean1d_sc(1,1,1,is), wrk1d, wrk2d, wrk3d)
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
     CALL DNS_SAVE_SCBDGIJ_P1(NNstat, m_p_x, m_p_y, m_p_z, &
          m_z1, z1(1,1,1,is), tmp7, tmp8, tmp9, &
          mean1d_sc(1,1,1,is), wrk2d)
  ENDDO

#undef m_p_x 
#undef m_p_y 
#undef m_p_z

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING DNS_SAVE_SCBDGIJ' )
#endif

  RETURN
END SUBROUTINE DNS_SAVE_SCBDGIJ
