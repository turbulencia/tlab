#include "dns_const.h"

!########################################################################
!# Skewsymmetric formulation according to Erlebacher, 1992.
!# Derived from RHS_FLOW_EULER_DIVERGENCE, 12 additional derivative
!# operations are added.
!# 27 derivative operations.
!########################################################################
subroutine RHS_FLOW_EULER_SKEWSYMMETRIC()
    use TLAB_CONSTANTS, only: wp, wi
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, inb_scal, imode_eqns
    use TLAB_VARS, only: g, buoyancy
    use TLAB_POINTERS
    use TLAB_ARRAYS, only: s
    use DNS_ARRAYS, only: hq, hs
    use THERMO_VARS, only: CRATIO_INV
    use OPR_PARTIAL

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i, is
    real(wp) g1, g2, g3, dummy

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

    bcs = 0

    g1 = buoyancy%vector(1)
    g2 = buoyancy%vector(2)
    g3 = buoyancy%vector(3)

! ###################################################################
! Terms \rho u in mass, u-momentum and energy equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*u(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i) + p(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp5)
    hq(:,5) = hq(:,5) - 2.0_wp*tmp5

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    do i = 1, imax*jmax*kmax
        hq(i,1) = hq(i,1) - (tmp2(i) + tmp3(i) + tmp4(i) + u(i)*tmp5(i)) + g1*rho(i)
        hq(i,2) = hq(i,2) - v(i)*tmp5(i)
        hq(i,3) = hq(i,3) - w(i)*tmp5(i)
    end do

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
    if (imode_eqns == DNS_EQNS_TOTAL) then
        hq(:,4) = hq(:,4) - (e + CRATIO_INV*0.5_wp*(u*u + v*v + w*w))*tmp5
! Internal energy
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        hq(:,4) = hq(:,4) - e*tmp5
    end if

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
    do is = 1, inb_scal
        hs(:, is) = hs(:, is) - s(:, is)*tmp5
    end do

! ###################################################################
! Terms \rho v in mass, v-momentum and energy equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*v(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i) + p(i)
        tmp1(i) = tmp4(i)*u(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp5)
    hq(:,5) = hq(:,5) - 2.0_wp*tmp5

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    do i = 1, imax*jmax*kmax
        hq(i,1) = hq(i,1) - u(i)*tmp5(i)
        hq(i,2) = hq(i,2) - (tmp2(i) + tmp3(i) + tmp4(i) + v(i)*tmp5(i)) + g2*rho(i)
        hq(i,3) = hq(i,3) - w(i)*tmp5(i)
    end do

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
    if (imode_eqns == DNS_EQNS_TOTAL) then
        hq(:,4) = hq(:,4) - (e + CRATIO_INV*0.5_wp*(u*u + v*v + w*w))*tmp5
! Internal energy
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        hq(:,4) = hq(:,4) - e*tmp5
    end if

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
    do is = 1, inb_scal
        hs(:, is) = hs(:, is) - s(:, is)*tmp5
    end do

! ###################################################################
! Terms \rho w in mass, w-momentum and energy equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*w(i)
        tmp3(i) = tmp4(i)*w(i) + p(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp5)
    hq(:,5) = hq(:,5) - 2.0_wp*tmp5

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    do i = 1, imax*jmax*kmax
        hq(i,1) = hq(i,1) - u(i)*tmp5(i)
        hq(i,2) = hq(i,2) - v(i)*tmp5(i)
        hq(i,3) = hq(i,3) - (tmp2(i) + tmp3(i) + tmp4(i) + w(i)*tmp5(i)) + g3*rho(i)
    end do

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
    if (imode_eqns == DNS_EQNS_TOTAL) then
        hq(:,4) = hq(:,4) - (e + CRATIO_INV*0.5_wp*(u*u + v*v + w*w))*tmp5
! Internal energy
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        hq(:,4) = hq(:,4) - e*tmp5
    end if

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
    do is = 1, inb_scal
        hs(:, is) = hs(:, is) - s(:, is)*tmp5
    end do

! ###################################################################
! Energy equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
    if (imode_eqns == DNS_EQNS_TOTAL) then
        do i = 1, imax*jmax*kmax
            dummy = 0.5_wp*rho(i)*(e(i) + CRATIO_INV*0.5_wp*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))) + CRATIO_INV*p(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
            tmp5(i) = (dummy - CRATIO_INV*p(i))/rho(i)
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
            dummy = 0.5_wp*rho(i)*e(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
            tmp5(i) = dummy/rho(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        hq(:,4) = hq(:,4) - (tmp2 + tmp3 + tmp4)
    end if

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! ###################################################################
! energy equation (array tmp5)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp5, tmp2)
    hq(:,4) = hq(:,4) - rho*(u*tmp2 + v*tmp3 + w*tmp4)

! u-momentum equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp2)
    hq(:,1) = hq(:,1) - 0.5_wp*rho*(u*tmp2 + v*tmp3 + w*tmp4)

! v-momentum equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp2)
    hq(:,2) = hq(:,2) - 0.5_wp*rho*(u*tmp2 + v*tmp3 + w*tmp4)

! w-momentum equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp2)
    hq(:,3) = hq(:,3) - 0.5_wp*rho*(u*tmp2 + v*tmp3 + w*tmp4)

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

    return
end subroutine RHS_FLOW_EULER_SKEWSYMMETRIC
