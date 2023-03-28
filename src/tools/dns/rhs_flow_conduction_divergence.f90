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
    use TLAB_CONSTANTS, only: wi, wp
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: g
    use TLAB_VARS, only: idiffusion, visc, prandtl
    use TLAB_POINTERS
    use TLAB_ARRAYS, only: s
    use DNS_ARRAYS, only: hq
    use THERMO_VARS, only: imixture
    use BOUNDARY_BCS
    use OPR_PARTIAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1)
    real(wp) cond

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

    bcs = 0

    if (idiffusion == EQNS_NONE) then
        cond = 0.0_wp
    else
        cond = visc/prandtl
    end if

#define diff_x(i) tmp5(i)
#define diff_y(i) tmp6(i)
#define diff_z(i) tmp7(i)

! calculate the enthalpy
    call THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, s, T, tmp4)

! enthalpy gradient
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)

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
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp1, tmp2)
    hq(:,4) = hq(:,4) + tmp2 + tmp3 + tmp4

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_CONDUCTION_DIVERGENCE')
#endif

    return
end subroutine RHS_FLOW_CONDUCTION_DIVERGENCE
