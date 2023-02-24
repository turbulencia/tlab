#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Computing terms depending on diffusion velocities, i.e. diffusion
!# term in the species eqns. and enthalpy transport in energy eqn.
!#
!# Using 2nd order derivative finite difference operators
!#
!########################################################################
subroutine RHS_SCAL_DIFFUSION_EXPLICIT(is, vis, z1, T, zh1, h4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)

    use TLAB_CONSTANTS, only: efile
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use TLAB_VARS, only: idiffusion, visc, prandtl, schmidt
    use THERMO_VARS, only: imixture, THERMO_AI, THERMO_TLIM, NSP, NCP
    use OPR_PARTIAL

    implicit none

    TINTEGER is
    TREAL, dimension(isize_field), intent(IN) :: vis, T
    TREAL, dimension(isize_field, *), intent(IN) :: z1
    TREAL, dimension(isize_field), intent(OUT) :: zh1, h4
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6

! -------------------------------------------------------------------
    TINTEGER bcs(2, 2), i, im, icp
    TREAL diff, cond

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

    bcs = 0

    if (idiffusion == EQNS_NONE) then; diff = C_0_R; cond = C_0_R
    else; diff = visc/schmidt(is); cond = visc/prandtl; end if

! ###################################################################
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), z1, tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), z1, tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), z1, tmp3, tmp4)
    zh1 = zh1 + diff*vis*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
    if (imixture > 0 .and. is < NSP .and. schmidt(is) /= prandtl) then
        tmp4 = C_0_R

        do i = 1, imax*jmax*kmax
            if (T(i) < THERMO_TLIM(3, is)) then; im = 2
            else; im = 1; end if
            do icp = NCP, 1, -1
                tmp4(i) = tmp4(i)*T(i) + (THERMO_AI(icp, im, is) - THERMO_AI(icp, im, NSP))/M_REAL(icp)
            end do
! factor (diff-cond) added now
            tmp4(i) = (diff - cond)*(tmp4(i)*T(i) + THERMO_AI(6, im, is) - THERMO_AI(6, im, NSP))
        end do
        h4 = h4 + vis*tmp4*(tmp1 + tmp2 + tmp3)

! cross-gradients
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), z1, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), z1, tmp5)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), z1, tmp6)
        h4 = h4 + vis*(tmp1*tmp4 + tmp2*tmp5 + tmp3*tmp6)

    end if

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_DIFFUSION_EXPLICIT')
#endif

    return
end subroutine RHS_SCAL_DIFFUSION_EXPLICIT
