#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

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
SUBROUTINE RHS_FLOW_CONDUCTION_EXPLICIT(vis, z1, T, h4, tmp1,tmp2,tmp3,tmp4,tmp5, wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : idiffusion, itransport, visc, prandtl
  USE BOUNDARY_BCS

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: vis, T
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: h4
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TREAL cond

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     CALL IO_WRITE_ASCII(efile,'RHS_FLOW_CONDUCTION_EXPLICIT. Only constant viscosity.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! ###################################################################
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; cond = C_0_R
  ELSE;                                  cond = visc/prandtl; ENDIF

! calculate the enthalpy
  CALL THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, z1, T, tmp4)

! total flux
  CALL OPR_PARTIAL_Z(OPR_P2, imax,jmax,kmax, bcs_out(1,1,3), g(3), tmp4, tmp3, tmp5, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2, imax,jmax,kmax, bcs_out(1,1,2), g(2), tmp4, tmp2, tmp5, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2, imax,jmax,kmax, bcs_out(1,1,1), g(1), tmp4, tmp1, tmp5, wrk2d,wrk3d)
  h4 = h4 + cond*vis*( tmp1 + tmp2 + tmp3 )

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_CONDUCTION_EXPLICIT
