#include "dns_const.h"

!########################################################################
!# 15 first derivative operations.
!########################################################################
subroutine RHS_FLOW_EULER_DIVERGENCE()
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, imode_eqns
    use TLAB_VARS, only: g, buoyancy
    use TLab_Pointers
    use DNS_ARRAYS, only: hq
    use Thermodynamics, only: CRATIO_INV
    use OPR_PARTIAL
    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i
    real(wp) g1, g2, g3, dummy

! ###################################################################
    bcs = 0

    g1 = buoyancy%vector(1)
    g2 = buoyancy%vector(2)
    g3 = buoyancy%vector(3)

! ###################################################################
! Terms \rho u in mass and u-momentum equations
! Explicit loops save memory calls
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = rho(i)*u(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i) + p(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp5)
    hq(:,5) = hq(:,5) - tmp5

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    hq(:,1) = hq(:,1) - (tmp2 + tmp3 + tmp4) + g1*rho

! ###################################################################
! Terms \rho v in mass and v-momentum equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = rho(i)*v(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i) + p(i)
        tmp1(i) = tmp4(i)*u(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp5)
    hq(:,5) = hq(:,5) - tmp5

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    hq(:,2) = hq(:,2) - (tmp2 + tmp3 + tmp4) + g2*rho

! ###################################################################
! Terms \rho w in mass and w-momentum equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = rho(i)*w(i)
        tmp3(i) = tmp4(i)*w(i) + p(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp5)
    hq(:,5) = hq(:,5) - tmp5

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    hq(:,3) = hq(:,3) - (tmp2 + tmp3 + tmp4) + g3*rho

! ###################################################################
! Total enery equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
    if (imode_eqns == DNS_EQNS_TOTAL) then
        do i = 1, imax*jmax*kmax
            dummy = rho(i)*(e(i) + CRATIO_INV*0.5_wp*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))) + CRATIO_INV*p(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        hq(:,4) = hq(:,4) - (tmp2 + tmp3 + tmp4) + CRATIO_INV*rho*(g1*u + g2*v + g3*w)

! -------------------------------------------------------------------
! Internal energy formulation
! the term p div u is add in RHS_FLOW_VISCOUS to save derivatives
! -------------------------------------------------------------------
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        do i = 1, imax*jmax*kmax
            dummy = rho(i)*e(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        hq(:,4) = hq(:,4) - (tmp2 + tmp3 + tmp4)

    end if

    return
end subroutine RHS_FLOW_EULER_DIVERGENCE
