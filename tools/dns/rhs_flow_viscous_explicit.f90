!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/30 - J.P. Mellado
!#              Created
!#
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
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE RHS_FLOW_VISCOUS_EXPLICIT(dx,dy,dz, vis, u,v,w,p, h1,h2,h3,h4, &
     tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL vis(*), u(*), v(*), w(*), p(*)
  TREAL h1(*), h2(*), h3(*), h4(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i1vsin, i1vsout, imxvsin, imxvsout
  TINTEGER j1vsin, j1vsout, jmxvsin, jmxvsout
  TINTEGER k1vsin, k1vsout, kmxvsin, kmxvsout
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

  prefactor = (gama0-C_1_R)*mach*mach
  c13 = C_1_R/C_3_R
  c23 = C_2_R/C_3_R

#include "dns_bcs_inf.h"
#include "dns_bcs_out.h"

! ###################################################################
! First derivatives in energy equation
! ###################################################################
  dummy = prefactor*visc
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + dummy*vis(i)*( (tmp2(i)+tmp3(i))*(tmp2(i)+tmp3(i)) )
  ENDDO
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + dummy*vis(i)*( (tmp2(i)+tmp3(i))*(tmp2(i)+tmp3(i)) )
  ENDDO
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + dummy*vis(i)*( (tmp2(i)+tmp3(i))*(tmp2(i)+tmp3(i)) )
  ENDDO

! ###################################################################
! Dilatation part
! ###################################################################
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)

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
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp1, i0,i0, i1vsin, imxvsin, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp2, i0,i0, j1vsout, jmxvsout, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp3, i0,i0, k1vsout, kmxvsout, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp5, tmp4, i1vsin, imxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h1(i) = h1(i) + vis(i)*visc*(tmp1(i) + tmp2(i) + tmp3(i) + c13*tmp4(i))
  ENDDO

  CALL PARTIAL_XX(i0, iunifx,imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp1, i0,i0, i1vsout, imxvsout, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp2, i0,i0, j1vsin, jmxvsin, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp3, i0,i0, k1vsout, kmxvsout, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp5, tmp4, j1vsin, jmxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h2(i) = h2(i) + vis(i)*visc*(tmp1(i) + tmp2(i) + tmp3(i) + c13*tmp4(i))
  ENDDO

  CALL PARTIAL_XX(i0, iunifx,imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp1, i0,i0, i1vsout, imxvsout, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp2, i0,i0, j1vsout, jmxvsout, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp3, i0,i0, k1vsin, kmxvsin, tmp4, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp5, tmp4, k1vsin, kmxvsin, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h3(i) = h3(i) + vis(i)*visc*(tmp1(i) + tmp2(i) + tmp3(i) + c13*tmp4(i))
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_VISCOUS_EXPLICIT
