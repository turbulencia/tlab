#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthapy transport in energy eqn. The
!# latter is stored cumulatively in arrays diff to be employed later
!# in the routine RHS_FLOW_CONDUCTION.
!#
!########################################################################
subroutine RHS_SCAL_DIFFUSION_DIVERGENCE &
    (is, vis, z1, T, zh1, diff_x, diff_y, diff_z, tmp1, tmp2, tmp3, tmp4)
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use TLAB_VARS, only: idiffusion, visc, prandtl, schmidt
    use THERMO_VARS, only: imixture, THERMO_AI, THERMO_TLIM, NSP, NCP
    use BOUNDARY_BCS
    use OPR_PARTIAL

    implicit none

    TINTEGER is
    TREAL, dimension(isize_field), intent(IN) :: vis, T
    TREAL, dimension(isize_field, *), intent(IN) :: z1
    TREAL, dimension(isize_field), intent(OUT) :: zh1
    TREAL, dimension(isize_field), intent(OUT) :: diff_x, diff_y, diff_z
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! -------------------------------------------------------------------
    TINTEGER bcs(2, 1)

    TINTEGER i, im, icp
    TREAL diff, cond, dummy
    TREAL ENTHALPY_L, ENTHALPY_G

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

    bcs = 0

    if (idiffusion == EQNS_NONE) then; diff = C_0_R; cond = C_0_R
    else; diff = visc/schmidt(is); cond = visc/prandtl; end if

! -------------------------------------------------------------------
! mass fraction gradients
! -------------------------------------------------------------------
! diffusion velocities in special cases
    if (imixture == MIXT_TYPE_AIRWATER) then
        tmp4 = (z1(:, 1) - z1(:, 2))/(C_1_R - z1(:, 2))
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)

! standard diffusion velocities
    else
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), z1, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), z1, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), z1, tmp3)

    end if

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
    if (imixture > 0 .and. is < NSP .and. schmidt(is) /= prandtl) then
        tmp4 = C_0_R

        do i = 1, imax*jmax*kmax
            if (T(i) < THERMO_TLIM(3, is)) then
                im = 2
            else
                im = 1
            end if
            do icp = NCP, 1, -1
                tmp4(i) = tmp4(i)*T(i) + (THERMO_AI(icp, im, is) - THERMO_AI(icp, im, NSP))/M_REAL(icp)
            end do
! factor (diff-cond) added now
            tmp4(i) = (diff - cond)*(tmp4(i)*T(i) + THERMO_AI(6, im, is) - THERMO_AI(6, im, NSP))
        end do
        diff_x = diff_x + tmp4*tmp1
        diff_y = diff_y + tmp4*tmp2
        diff_z = diff_z + tmp4*tmp3

    end if

! -------------------------------------------------------------------
! mass fraction transport by diffusion velocities
! -------------------------------------------------------------------
    tmp1 = diff*vis*tmp1
    tmp2 = diff*vis*tmp2
    tmp3 = diff*vis*tmp3

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp1, tmp2)
    zh1 = zh1 + tmp2 + tmp3 + tmp4

! ###################################################################
! Numerical liquid correction term in the AIRWATER case
! ###################################################################
    if (imixture == MIXT_TYPE_AIRWATER) then
        if (idiffusion == EQNS_NONE) then; diff = C_0_R
        else; diff = visc/schmidt(3); end if

! gradient of liquid content
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), z1(1, 2), tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), z1(1, 2), tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), z1(1, 2), tmp3)

! enthalpy equation
        do i = 1, imax*jmax*kmax
            ENTHALPY_G = &
                ((C_1_R - z1(i, 1))*(THERMO_AI(6, 1, 2) + THERMO_AI(1, 1, 2)*T(i)) &
                 + (z1(i, 1) - z1(i, 2))*(THERMO_AI(6, 1, 1) + THERMO_AI(1, 1, 1)*T(i)) &
                 )/(C_1_R - z1(i, 2))
            ENTHALPY_L = &
                THERMO_AI(6, 1, 3) + THERMO_AI(1, 1, 3)*T(i)
            diff_x(i) = diff_x(i) + diff*(ENTHALPY_L - ENTHALPY_G)*tmp1(i)
            diff_y(i) = diff_y(i) + diff*(ENTHALPY_L - ENTHALPY_G)*tmp2(i)
            diff_z(i) = diff_z(i) + diff*(ENTHALPY_L - ENTHALPY_G)*tmp3(i)
        end do

! scalar equation
        do i = 1, imax*jmax*kmax
            dummy = (C_1_R - z1(i, 1))/(C_1_R - z1(i, 2))
            tmp1(i) = diff*vis(i)*tmp1(i)*dummy
            tmp2(i) = diff*vis(i)*tmp2(i)*dummy
            tmp3(i) = diff*vis(i)*tmp3(i)*dummy
        end do

        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp1, tmp2)
        zh1 = zh1 + tmp2 + tmp3 + tmp4

    end if

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_DIVERGENCE')
#endif

    return
end subroutine RHS_SCAL_DIFFUSION_DIVERGENCE
