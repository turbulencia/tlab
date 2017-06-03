#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate div(tau) in terms of second order finite differences,
!# assuming constant viscosity.
!# The dissipation function is implemented only for the case of the
!# internal energy formulation.
!# The BCs are not exactly those imposed in the divergence formulation
!# because of the rearrangement into the dilatation part of part of
!# the second derivative. No strong impact has been observed due to this.
!# 12 first derivative operations and 9 second derivative operations.
!# 
!########################################################################
SUBROUTINE RHS_FLOW_VISCOUS_EXPLICIT(vis, u,v,w,p, h1,h2,h3,h4, tmp1,tmp2,tmp3,tmp4,tmp5, wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field, imode_eqns
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : itransport, visc, mach
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL,     ONLY : bcs_out, bcs_inf

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field), INTENT(IN)    :: vis, u,v,w,p
  TREAL, DIMENSION(isize_field), INTENT(OUT)   :: h1,h2,h3,h4
  TREAL, DIMENSION(isize_field), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),           INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1)
  TINTEGER i
  TREAL prefactor, c13, c23, dummy

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_VISCOUS_EXPLICIT. Only constant viscosity.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_VISCOUS_EXPLICIT. No total energy formulation.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  bcs = 0
  
  prefactor = (gama0-C_1_R)*mach*mach
  c13 = C_1_R/C_3_R
  c23 = C_2_R/C_3_R

! ###################################################################
! First derivatives in energy equation
! ###################################################################
  dummy = prefactor*visc
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp3, wrk3d, wrk2d,wrk3d)
  h4 = h4 + dummy*vis*( (tmp2+tmp3)*(tmp2+tmp3) )

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp3, wrk3d, wrk2d,wrk3d)
  h4 = h4 + dummy*vis*( (tmp2+tmp3)*(tmp2+tmp3) )

  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, tmp3, wrk3d, wrk2d,wrk3d)
  h4 = h4 + dummy*vis*( (tmp2+tmp3)*(tmp2+tmp3) )

! ###################################################################
! Dilatation part
! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
  dummy = c23*visc
  DO i = 1,imax*jmax*kmax
     tmp5(i) = tmp1(i) + tmp2(i) + tmp3(i)
     h4(i) = h4(i) + prefactor*( dummy*vis(i)*( &
          (tmp1(i)-tmp2(i))*(tmp1(i)-tmp2(i)) + &
          (tmp2(i)-tmp3(i))*(tmp2(i)-tmp3(i)) + &
          (tmp3(i)-tmp1(i))*(tmp3(i)-tmp1(i)) ) - p(i)*tmp5(i) )
  ENDDO

! ###################################################################
! Laplacian terms in the momentum equation
! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs_inf(1,1,1), g(1), u,    tmp1, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs_out(1,1,2), g(2), u,    tmp2, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs_out(1,1,3), g(3), u,    tmp3, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,1), g(1), tmp5, tmp4, wrk3d, wrk2d,wrk3d)
  h1 = h1 + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs_out(1,1,1), g(1), v,    tmp1, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs_inf(1,1,2), g(2), v,    tmp2, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs_out(1,1,3), g(3), v,    tmp3, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,2), g(2), tmp5, tmp4, wrk3d, wrk2d,wrk3d)
  h2 = h2 + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs_out(1,1,1), g(1), w,    tmp1, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs_out(1,1,2), g(2), w,    tmp2, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs_inf(1,1,3), g(3), w,    tmp3, tmp4,  wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_inf(1,2,3), g(3), tmp5, tmp4, wrk3d, wrk2d,wrk3d)
  h3 = h3 + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_VISCOUS_EXPLICIT
