#include "dns_const.h"

!########################################################################
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthalpy transport in energy eqn.
!#
!# Using 2nd order derivative finite difference operators
!########################################################################
subroutine RHS_SCAL_DIFFUSION_EXPLICIT(is)
    use TLab_Constants, only: efile, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
    use TLab_Memory, only: imax, jmax, kmax
    use FDM, only: g
    use NavierStokes, only: nse_diffusion
    use NavierStokes, only: visc, prandtl, schmidt
    use Thermodynamics, only: imixture, THERMO_AI, THERMO_TLIM, NSP, NCP
    use TLab_Pointers
    use TLab_Arrays, only: s
    use DNS_ARRAYS, only: hs, hq
    use OPR_PARTIAL
    implicit none

    integer, intent(in) :: is

! -------------------------------------------------------------------
    integer(wi) bcs(2, 2), i, im, icp
    real(wp) diff, cond

! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

    bcs = 0

    if (nse_diffusion == EQNS_NONE) then
        diff = 0.0_wp; cond = 0.0_wp
    else
        diff = visc/schmidt(is); cond = visc/prandtl
    end if

! ###################################################################
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), s(:, is), tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), s(:, is), tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), s(:, is), tmp3, tmp4)
    hs(:, is) = hs(:, is) + diff*vis*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
    if (imixture > 0 .and. is < NSP .and. schmidt(is) /= prandtl) then
        tmp4 = 0.0_wp

        do i = 1, imax*jmax*kmax
            if (T(i) < THERMO_TLIM(3, is)) then; im = 2
            else; im = 1; end if
            do icp = NCP, 1, -1
                tmp4(i) = tmp4(i)*T(i) + (THERMO_AI(icp, im, is) - THERMO_AI(icp, im, NSP))/real(icp, wp)
            end do
! factor (diff-cond) added now
            tmp4(i) = (diff - cond)*(tmp4(i)*T(i) + THERMO_AI(6, im, is) - THERMO_AI(6, im, NSP))
        end do
        hq(:, 4) = hq(:, 4) + vis*tmp4*(tmp1 + tmp2 + tmp3)

! cross-gradients
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(:, is), tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), tmp5)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(:, is), tmp6)
        hq(:, 4) = hq(:, 4) + vis*(tmp1*tmp4 + tmp2*tmp5 + tmp3*tmp6)

    end if

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

    return
end subroutine RHS_SCAL_DIFFUSION_EXPLICIT
