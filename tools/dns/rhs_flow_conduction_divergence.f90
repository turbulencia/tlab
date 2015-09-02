#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/06/08 - J.P. Mellado
!#              Created
!#
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
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE RHS_FLOW_CONDUCTION_DIVERGENCE&
     (dx, dy, dz, vis, z1, T, h4, diff_x, diff_y, diff_z,&
     tmp1, tmp2, tmp3, tmp4, wrk1d, wrk2d, wrk3d)

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)

  TREAL vis(*), T(*), z1(imax*jmax*kmax,*)
  TREAL diff_x(*), diff_y(*), diff_z(*)
  TREAL h4(*)
  TREAL tmp1(*), tmp2(*), tmp3(*), tmp4(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER i1vsout, imxvsout
  TINTEGER j1vsout, jmxvsout
  TINTEGER k1vsout, kmxvsout
  TINTEGER i
  TREAL cond

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

#include "dns_bcs_out.h"

! ###################################################################
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; cond = C_0_R
  ELSE;                                  cond = visc/prandtl; ENDIF

! calculate the enthalpy
  CALL THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, z1, T, tmp4)

! enthalpy gradient
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp4, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp4, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp4, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)

! Add diffusion velocity terms as required
  IF ( imixture .GT. 0 ) THEN
     DO i = 1,imax*jmax*kmax
        tmp1(i) = vis(i)*(cond*tmp1(i) + diff_x(i))
        tmp2(i) = vis(i)*(cond*tmp2(i) + diff_y(i))
        tmp3(i) = vis(i)*(cond*tmp3(i) + diff_z(i))
     ENDDO
  ELSE
     DO i = 1,imax*jmax*kmax
        tmp1(i) = cond*vis(i)*tmp1(i)
        tmp2(i) = cond*vis(i)*tmp2(i)
        tmp3(i) = cond*vis(i)*tmp3(i)
     ENDDO
  ENDIF

! total flux
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, tmp3, tmp4, k1vsout, kmxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, tmp2, tmp3, j1vsout, jmxvsout, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, tmp1, tmp2, i1vsout, imxvsout, wrk1d, wrk2d, wrk3d)
  DO i = 1,imax*jmax*kmax
     h4(i) = h4(i) + tmp2(i) + tmp3(i) + tmp4(i)
  ENDDO

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

  RETURN
END SUBROUTINE RHS_FLOW_CONDUCTION_DIVERGENCE
