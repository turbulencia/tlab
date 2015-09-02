!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/21 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Compute heat flux term in the energy equation
!#
!# div ( \mu/Pr grad h )
!#
!# using 2nd order derivative finite difference operator.
!# Mass diffusion contribbution in RHS_SCAL_DIFFSUION_EXPLICIT.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE RHS_FLOW_CONDUCTION_EXPLICIT&
     (dx,dy,dz, vis, z1, T, h4, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL vis(*), T(*), z1(imax*jmax*kmax,*)
  TREAL h4(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*), tmp5(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i1vsout, imxvsout
  TINTEGER j1vsout, jmxvsout
  TINTEGER k1vsout, kmxvsout
  TINTEGER i
  TREAL cond

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_CONDUCTION_EXPLICIT. Only constant viscosity.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

#include "dns_bcs_out.h"

! ###################################################################
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; cond = C_0_R
  ELSE;                                  cond = visc/prandtl; ENDIF

! calculate the enthalpy
  CALL THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, z1, T, tmp4)

! total flux
  CALL PARTIAL_ZZ(i0, iunifz, imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp4, tmp3, i0,i0, k1vsout, kmxvsout, tmp5, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_YY(i0, iunify, imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp4, tmp2, i0,i0, j1vsout, jmxvsout, tmp5, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_XX(i0, iunifx, imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp1, i0,i0, i1vsout, imxvsout, tmp5, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + cond*vis(i)*( tmp1(i) + tmp2(i) + tmp3(i) )
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_CONDUCTION_EXPLICIT
