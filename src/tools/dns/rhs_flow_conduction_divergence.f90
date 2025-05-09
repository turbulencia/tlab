#include "dns_const.h"

!########################################################################
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
subroutine RHS_FLOW_CONDUCTION_DIVERGENCE()
    use TLab_Constants, only: wi, wp
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
    use TLab_Memory, only: imax, jmax, kmax
    use FDM, only: g
    use NavierStokes, only: nse_diffusion, EQNS_NONE
    use NavierStokes, only: visc, prandtl
    use TLab_Pointers
    use TLab_Arrays, only: s
    use Thermodynamics, only: imixture
    use THERMO_CALORIC
    use DNS_ARRAYS, only: hq
    use BOUNDARY_BCS
    use OPR_Partial

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1)
    real(wp) cond

! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

    bcs = 0

    if (nse_diffusion == EQNS_NONE) then
        cond = 0.0_wp
    else
        cond = visc/prandtl
    end if

#define diff_x(i) tmp5(i)
#define diff_y(i) tmp6(i)
#define diff_z(i) tmp7(i)

! calculate the enthalpy
    call THERMO_CALORIC_ENTHALPY(imax*jmax*kmax, s, T, tmp4)

! enthalpy gradient
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)

! Add diffusion velocity terms as required
    if (imixture > 0) then
        tmp1 = vis*(cond*tmp1 + diff_x(:))
        tmp2 = vis*(cond*tmp2 + diff_y(:))
        tmp3 = vis*(cond*tmp3 + diff_z(:))

    else
        tmp1 = cond*vis*tmp1
        tmp2 = cond*vis*tmp2
        tmp3 = cond*vis*tmp3

    end if

! total flux
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp3, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp2, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp1, tmp2)
    hq(:,4) = hq(:,4) + tmp2 + tmp3 + tmp4

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

    return
end subroutine RHS_FLOW_CONDUCTION_DIVERGENCE
