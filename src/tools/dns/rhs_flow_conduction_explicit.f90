#include "dns_const.h"

!########################################################################
!# Compute heat flux term in the energy equation
!#
!# div ( \mu/Pr grad h )
!#
!# using 2nd order derivative finite difference operator.
!# Mass diffusion contribbution in RHS_SCAL_DIFFSUION_EXPLICIT.
!########################################################################
subroutine RHS_FLOW_CONDUCTION_EXPLICIT()
    use TLab_Constants, only: efile, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax
    use FDM, only: g
    use TLAB_VARS, only: idiffusion, visc, prandtl
    use TLab_Pointers
    use TLab_Arrays, only: s
    use THERMO_CALORIC
    use DNS_ARRAYS, only: hq
    use BOUNDARY_BCS
    use OPR_PARTIAL
    implicit none

! -------------------------------------------------------------------
    real(wp) cond

! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

! ###################################################################
    if (idiffusion == EQNS_NONE) then
        cond = 0.0_wp
    else
        cond = visc/prandtl
    end if

! calculate the enthalpy
    call THERMO_CALORIC_ENTHALPY(imax*jmax*kmax, s, T, tmp4)

! total flux
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp4, tmp3, tmp5)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp4, tmp2, tmp5)
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp4, tmp1, tmp5)
    hq(:,4) = hq(:,4) + cond*vis*(tmp1 + tmp2 + tmp3)

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_EXPLICIT')
#endif

    return
end subroutine RHS_FLOW_CONDUCTION_EXPLICIT
