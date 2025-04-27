#include "dns_const.h"

!########################################################################
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthapy transport in energy eqn. The
!# latter is stored cumulatively in arrays diff to be employed later
!# in the routine RHS_FLOW_CONDUCTION.
!########################################################################
subroutine RHS_SCAL_DIFFUSION_DIVERGENCE(is)
    use TLab_Constants, only: wi, wp
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
    use TLab_Memory, only: imax, jmax, kmax
    use FDM, only: g
    use NavierStokes, only: nse_diffusion
    use NavierStokes, only: visc, prandtl, schmidt
    use TLab_Pointers
    use TLab_Arrays, only: s
    use DNS_ARRAYS, only: hs
    use Thermodynamics, only: imixture, MIXT_TYPE_AIRWATER, THERMO_AI, THERMO_TLIM, NSP, NCP
    use OPR_Partial
    use BOUNDARY_BCS
    implicit none

    integer, intent(in) :: is

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i, im, icp
    real(wp) diff, cond, dummy
    real(wp) ENTHALPY_L, ENTHALPY_G

! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

    bcs = 0

    if (nse_diffusion == EQNS_NONE) then
        diff = 0.0_wp; cond = 0.0_wp
    else
        diff = visc/schmidt(is); cond = visc/prandtl
    end if

#define diff_x(i) tmp5(i)
#define diff_y(i) tmp6(i)
#define diff_z(i) tmp7(i)

! -------------------------------------------------------------------
! mass fraction gradients
! -------------------------------------------------------------------
! diffusion velocities in special cases
    if (imixture == MIXT_TYPE_AIRWATER) then
        tmp4 = (s(:, 1) - s(:, 2))/(1.0_wp - s(:, 2))
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)

! standard diffusion velocities
    else
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(:, is), tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(:, is), tmp3)

    end if

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
    if (imixture > 0 .and. is < NSP .and. schmidt(is) /= prandtl) then
        tmp4 = 0.0_wp

        do i = 1, imax*jmax*kmax
            if (T(i) < THERMO_TLIM(3, is)) then
                im = 2
            else
                im = 1
            end if
            do icp = NCP, 1, -1
                tmp4(i) = tmp4(i)*T(i) + (THERMO_AI(icp, im, is) - THERMO_AI(icp, im, NSP))/real(icp, wp)
            end do
! factor (diff-cond) added now
            tmp4(i) = (diff - cond)*(tmp4(i)*T(i) + THERMO_AI(6, im, is) - THERMO_AI(6, im, NSP))
        end do
        diff_x(:) = diff_x(:) + tmp4*tmp1
        diff_y(:) = diff_y(:) + tmp4*tmp2
        diff_z(:) = diff_z(:) + tmp4*tmp3

    end if

! -------------------------------------------------------------------
! mass fraction transport by diffusion velocities
! -------------------------------------------------------------------
    tmp1 = diff*vis*tmp1
    tmp2 = diff*vis*tmp2
    tmp3 = diff*vis*tmp3

    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp3, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp2, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp1, tmp2)
    hs(:, is) = hs(:, is) + tmp2 + tmp3 + tmp4

! ###################################################################
! Numerical liquid correction term in the AIRWATER case
! ###################################################################
    if (imixture == MIXT_TYPE_AIRWATER) then
        if (nse_diffusion == EQNS_NONE) then; diff = 0.0_wp
        else; diff = visc/schmidt(3); end if

! gradient of liquid content
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, 2), tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, 2), tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, 2), tmp3)

! enthalpy equation
        do i = 1, imax*jmax*kmax
            ENTHALPY_G = &
                ((1.0_wp - s(i, 1))*(THERMO_AI(6, 1, 2) + THERMO_AI(1, 1, 2)*T(i)) &
                 + (s(i, 1) - s(i, 2))*(THERMO_AI(6, 1, 1) + THERMO_AI(1, 1, 1)*T(i)) &
                 )/(1.0_wp - s(i, 2))
            ENTHALPY_L = &
                THERMO_AI(6, 1, 3) + THERMO_AI(1, 1, 3)*T(i)
            diff_x(i) = diff_x(i) + diff*(ENTHALPY_L - ENTHALPY_G)*tmp1(i)
            diff_y(i) = diff_y(i) + diff*(ENTHALPY_L - ENTHALPY_G)*tmp2(i)
            diff_z(i) = diff_z(i) + diff*(ENTHALPY_L - ENTHALPY_G)*tmp3(i)
        end do

! scalar equation
        do i = 1, imax*jmax*kmax
            dummy = (1.0_wp - s(i, 1))/(1.0_wp - s(i, 2))
            tmp1(i) = diff*vis(i)*tmp1(i)*dummy
            tmp2(i) = diff*vis(i)*tmp2(i)*dummy
            tmp3(i) = diff*vis(i)*tmp3(i)*dummy
        end do

        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp3, tmp4)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp2, tmp3)
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp1, tmp2)
        hs(:, is) = hs(:, is) + tmp2 + tmp3 + tmp4

    end if

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

    return
end subroutine RHS_SCAL_DIFFUSION_DIVERGENCE
