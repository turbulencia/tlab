#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate div(tau) in terms of second order finite differences,
!# assuming constant viscosity.
!# The dissipation function is implemented only for the case of the
!# internal energy formulation.
!# The BCs are not exactly those imposed in the divergence formulation
!# because of the rearrangement into the dilatation part of part of
!# the second derivative. No strong impact has been observed due to this.
!# 12 first derivative operations and 9 second derivative operations.
!#
!########################################################################
subroutine RHS_FLOW_VISCOUS_EXPLICIT(vis, u, v, w, p, h1, h2, h3, h4, tmp1, tmp2, tmp3, tmp4, tmp5)

    use TLAB_CONSTANTS, only: efile
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use TLAB_VARS, only: visc, mach
    use THERMO_VARS, only: gama0
    use BOUNDARY_BCS
    use OPR_PARTIAL

    implicit none

    TREAL, dimension(isize_field), intent(IN) :: vis, u, v, w, p
    TREAL, dimension(isize_field), intent(OUT) :: h1, h2, h3, h4
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5

! -------------------------------------------------------------------
    TINTEGER bcs(2, 1)
    TINTEGER i
    TREAL prefactor, c13, c23, dummy

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

    bcs = 0

    prefactor = (gama0 - C_1_R)*mach*mach
    c13 = C_1_R/C_3_R
    c23 = C_2_R/C_3_R

! ###################################################################
! First derivatives in energy equation
! ###################################################################
    dummy = prefactor*visc
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp2)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp3)
    h4 = h4 + dummy*vis*((tmp2 + tmp3)*(tmp2 + tmp3))

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp2)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp3)
    h4 = h4 + dummy*vis*((tmp2 + tmp3)*(tmp2 + tmp3))

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp2)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp3)
    h4 = h4 + dummy*vis*((tmp2 + tmp3)*(tmp2 + tmp3))

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
        h4(i) = h4(i) + prefactor*(dummy*vis(i)*( &
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
    h1 = h1 + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), v, tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), v, tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), v, tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), tmp5, tmp4)
    h2 = h2 + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), w, tmp1, tmp4)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), w, tmp2, tmp4)
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), w, tmp3, tmp4)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), tmp5, tmp4)
    h3 = h3 + vis*visc*(tmp1 + tmp2 + tmp3 + c13*tmp4)

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_EXPLICIT')
#endif

    return
end subroutine RHS_FLOW_VISCOUS_EXPLICIT
