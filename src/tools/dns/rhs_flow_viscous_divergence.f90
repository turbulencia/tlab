#include "dns_const.h"

!########################################################################
!# The case of internal energy adds the term p div u here to avoid the
!# computation of the dilatation in RHS_FLOW_EULER_?, so it is
!# not only the viscous part in that case...
!# Internal energy eqn formulation does 18 derivatives.
!# Total energy eqn formulation does 21 derivatives.
!########################################################################
subroutine RHS_FLOW_VISCOUS_DIVERGENCE()
    use TLAB_CONSTANTS, only: wp, wi
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLab_WorkFlow, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, imode_eqns
    use TLAB_VARS, only: g
    use TLAB_VARS, only: visc
    use TLAB_POINTERS
    use DNS_ARRAYS, only: hq
    use Thermodynamics, only: CRATIO_INV
    use BOUNDARY_BCS
    use OPR_PARTIAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i
    real(wp) dil, c23

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_VISCOUS_DIVERGENCE')
#endif

    bcs = 0
    c23 = 2.0_wp/3.0_wp*visc

#define tau_xx(i) tmp4(i)
#define tau_xy(i) tmp5(i)
#define tau_xz(i) tmp6(i)
#define tau_yy(i) tmp7(i)
#define tau_yz(i) tmp8(i)
#define tau_zz(i) tmp9(i)

! ###################################################################
! Define viscous stress tensor
! Add corresponding terms in the internal energy equation if needed
! ###################################################################
! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)
    if (imode_eqns == DNS_EQNS_INTERNAL) then ! internal energy equation
        do i = 1, imax*jmax*kmax
            tau_xx(i) = c23*vis(i)*(2.0_wp*tmp1(i) - (tmp2(i) + tmp3(i)))
            tau_yy(i) = c23*vis(i)*(2.0_wp*tmp2(i) - (tmp1(i) + tmp3(i)))
            tau_zz(i) = c23*vis(i)*(2.0_wp*tmp3(i) - (tmp1(i) + tmp2(i)))
            dil = tmp1(i) + tmp2(i) + tmp3(i)
            hq(i, 4) = hq(i, 4) + CRATIO_INV*(tau_xx(i)*tmp1(i) + tau_yy(i)*tmp2(i) + tau_zz(i)*tmp3(i) - p(i)*dil)
        end do

    else if (imode_eqns == DNS_EQNS_TOTAL) then ! total energy equation
        do i = 1, imax*jmax*kmax
            tau_xx(i) = c23*vis(i)*(2.0_wp*tmp1(i) - (tmp2(i) + tmp3(i)))
            tau_yy(i) = c23*vis(i)*(2.0_wp*tmp2(i) - (tmp1(i) + tmp3(i)))
            tau_zz(i) = c23*vis(i)*(2.0_wp*tmp3(i) - (tmp1(i) + tmp2(i)))
        end do

    end if

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp2)
    if (imode_eqns == DNS_EQNS_INTERNAL) then ! internal energy equation
        do i = 1, imax*jmax*kmax
            tau_xy(i) = visc*vis(i)*(tmp1(i) + tmp2(i))
            hq(i, 4) = hq(i, 4) + CRATIO_INV*(tau_xy(i)*(tmp1(i) + tmp2(i)))
        end do

    else if (imode_eqns == DNS_EQNS_TOTAL) then ! total energy equation
        tau_xy(:) = visc*vis*(tmp1 + tmp2)

    end if

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp1)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp2)

    if (imode_eqns == DNS_EQNS_INTERNAL) then ! internal energy equation
        do i = 1, imax*jmax*kmax
            tau_xz(i) = visc*vis(i)*(tmp1(i) + tmp2(i))
            hq(i, 4) = hq(i, 4) + CRATIO_INV*(tau_xz(i)*(tmp1(i) + tmp2(i)))
        end do

    else if (imode_eqns == DNS_EQNS_TOTAL) then ! total energy equation
        tau_xz(:) = visc*vis*(tmp1 + tmp2)

    end if

    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp1)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp2)

    if (imode_eqns == DNS_EQNS_INTERNAL) then ! internal energy equation
        do i = 1, imax*jmax*kmax
            tau_yz(i) = visc*vis(i)*(tmp1(i) + tmp2(i))
            hq(i, 4) = hq(i, 4) + CRATIO_INV*(tau_yz(i)*(tmp1(i) + tmp2(i)))
        end do

    else if (imode_eqns == DNS_EQNS_TOTAL) then ! total energy equation
        tau_yz(:) = visc*vis*(tmp1 + tmp2)

    end if

! ###################################################################
! Momentum equation
! ###################################################################
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), tau_xx(:), tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tau_xy(:), tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tau_xz(:), tmp3)
    hq(:, 1) = hq(:, 1) + tmp1 + tmp2 + tmp3

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tau_xy(:), tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), tau_yy(:), tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tau_yz(:), tmp3)
    hq(:, 2) = hq(:, 2) + tmp1 + tmp2 + tmp3

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tau_xz(:), tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tau_yz(:), tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), tau_zz(:), tmp3)
    hq(:, 3) = hq(:, 3) + tmp1 + tmp2 + tmp3

! ###################################################################
! Energy equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
    if (imode_eqns == DNS_EQNS_TOTAL) then
        do i = 1, imax*jmax*kmax
            tau_xx(i) = tau_xx(i)*u(i) + tau_xy(i)*v(i) + tau_xz(i)*w(i)
            tau_yy(i) = tau_xy(i)*u(i) + tau_yy(i)*v(i) + tau_yz(i)*w(i)
            tau_zz(i) = tau_xz(i)*u(i) + tau_yz(i)*v(i) + tau_zz(i)*w(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tau_xx(:), tmp3)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tau_yy(:), tmp2)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tau_zz(:), tmp1)
        hq(:, 4) = hq(:, 4) + CRATIO_INV*(tmp1 + tmp2 + tmp3)

    end if

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_VISCOUS_DIVERGENCE')
#endif

    return
end subroutine RHS_FLOW_VISCOUS_DIVERGENCE
