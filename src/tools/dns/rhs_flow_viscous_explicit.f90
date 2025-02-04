#include "dns_const.h"

!########################################################################
!# Calculate div(tau) in terms of second order finite differences,
!# assuming constant viscosity.
!# The dissipation function is implemented only for the case of the
!# internal energy formulation.
!# The BCs are not exactly those imposed in the divergence formulation
!# because of the rearrangement into the dilatation part of part of
!# the second derivative. No strong impact has been observed due to this.
!# 12 first derivative operations and 9 second derivative operations.
!########################################################################
subroutine RHS_FLOW_VISCOUS_EXPLICIT()
    use TLab_Constants, only: efile, wi, wp
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
    use TLab_Memory, only: imax, jmax, kmax
    use FDM, only: g
    use NavierStokes, only: visc
    use TLab_Pointers
    use DNS_ARRAYS, only: hq
    use Thermodynamics, only: CRATIO_INV
    use BOUNDARY_BCS
    use OPR_PARTIAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i
    real(wp) c13, c23, dummy

! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

    bcs = 0

    c13 = 1.0_wp/3.0_wp
    c23 = 2.0_wp/3.0_wp

! ###################################################################
! First derivatives in energy equation
! ###################################################################
    dummy = CRATIO_INV*visc
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp2)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp3)
    hq(:, 4) = hq(:, 4) + dummy*vis*((tmp2 + tmp3)*(tmp2 + tmp3))

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp2)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp3)
    hq(:, 4) = hq(:, 4) + dummy*vis*((tmp2 + tmp3)*(tmp2 + tmp3))

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp2)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp3)
    hq(:, 4) = hq(:, 4) + dummy*vis*((tmp2 + tmp3)*(tmp2 + tmp3))

! ###################################################################
! Dilatation part
! ###################################################################
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
    dummy = c23*visc
    do i = 1, imax*jmax*kmax
        tmp5(i) = tmp1(i) + tmp2(i) + tmp3(i)
        hq(i, 4) = hq(i, 4) + CRATIO_INV*(dummy*vis(i)*( &
                                         (tmp1(i) - tmp2(i))*(tmp1(i) - tmp2(i)) + &
                                         (tmp2(i) - tmp3(i))*(tmp2(i) - tmp3(i)) + &
                                         (tmp3(i) - tmp1(i))*(tmp3(i) - tmp1(i))) - p(i)*tmp5(i))
    end do

! ###################################################################
! Laplacian terms in the momentum equation
! ###################################################################
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), u, tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), u, tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), u, tmp3, tmp4)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), tmp5, tmp4)
    hq(:, 1) = hq(:, 1) + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), v, tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), v, tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), v, tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), tmp5, tmp4)
    hq(:, 2) = hq(:, 2) + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), w, tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), w, tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), w, tmp3, tmp4)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), tmp5, tmp4)
    hq(:, 3) = hq(:, 3) + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

    return
end subroutine RHS_FLOW_VISCOUS_EXPLICIT
