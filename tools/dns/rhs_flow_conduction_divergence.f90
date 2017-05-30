#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Compute heat flux term in the energy equation
!# div ( \mu/Pr grad h + \sum \mu(1/Sc_i-1/Pr)(h_i-h_N) grad Y_i )
!#
!# The contribution from enthalpy transport by diffusion velocities 
!# enters through arrays diff_i. 
!# Obviously, RHS_SCAL_DIFFUSION must precede this routine, so that 
!# diff arrays are filled.
!#
!########################################################################
SUBROUTINE RHS_FLOW_CONDUCTION_DIVERGENCE&
     (vis, z1, T, h4, diff_x,diff_y,diff_z, tmp1,tmp2,tmp3,tmp4, wrk2d,wrk3d)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : idiffusion, visc, prandtl
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL,     ONLY : bcs_out

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: vis, T
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: z1
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: h4
  TREAL, DIMENSION(isize_field),   INTENT(IN)    :: diff_x, diff_y, diff_z
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,1)
  TREAL cond

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

  bcs = 0
  
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; cond = C_0_R
  ELSE;                                  cond = visc/prandtl; ENDIF

! calculate the enthalpy
  CALL THERMO_CALORIC_ENTHALPY(imax,jmax,kmax, z1, T, tmp4)

! enthalpy gradient
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp4, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d,wrk3d)

! Add diffusion velocity terms as required
  IF ( imixture .GT. 0 ) THEN
     tmp1 = vis *( cond *tmp1 +diff_x )
     tmp2 = vis *( cond *tmp2 +diff_y )
     tmp3 = vis *( cond *tmp3 +diff_z )
     
  ELSE
     tmp1 = cond *vis *tmp1
     tmp2 = cond *vis *tmp2
     tmp3 = cond *vis *tmp3
     
  ENDIF

! total flux
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs_out(1,2,3), g(3), tmp3, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_out(1,2,2), g(2), tmp2, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs_out(1,2,1), g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  h4 = h4 +tmp2 +tmp3 +tmp4

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_CONDUCTION_DIVERGENCE
