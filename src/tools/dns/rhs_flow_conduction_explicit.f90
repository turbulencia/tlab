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
subroutine RHS_FLOW_CONDUCTION_EXPLICIT(vis, z1, T, h4, tmp1, tmp2, tmp3, tmp4, tmp5)

    use TLAB_CONSTANTS, only: efile
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use TLAB_VARS, only: idiffusion, visc, prandtl
    use BOUNDARY_BCS
    use OPR_PARTIAL
    implicit none

    TREAL, dimension(isize_field), intent(IN) :: vis, T
    TREAL, dimension(isize_field, *), intent(IN) :: z1
    TREAL, dimension(isize_field), intent(OUT) :: h4
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5

! -------------------------------------------------------------------
    TREAL cond

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

! ###################################################################
    if (idiffusion == EQNS_NONE) then; cond = C_0_R
    else; cond = visc/prandtl; end if

! calculate the enthalpy
    call THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, z1, T, tmp4)

! total flux
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp4, tmp3, tmp5)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp4, tmp2, tmp5)
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp4, tmp1, tmp5)
    h4 = h4 + cond*vis*(tmp1 + tmp2 + tmp3)

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

    return
end subroutine RHS_FLOW_CONDUCTION_EXPLICIT
