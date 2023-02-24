#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Skewsymmetric formulation according to Erlebacher, 1992.
!# Derived from RHS_FLOW_EULER_DIVERGENCE, 12 additional derivative
!# operations are added.
!# 27 derivative operations.
!#
!########################################################################
subroutine RHS_FLOW_EULER_SKEWSYMMETRIC(rho, u, v, w, p, e, z1, h0, h1, h2, h3, h4, zh1, &
                                        tmp1, tmp2, tmp3, tmp4, tmp5)
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif

    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, imode_eqns
    use TLAB_VARS, only: g, buoyancy
    use TLAB_VARS, only: mach
    use THERMO_VARS, only: gama0
    use OPR_PARTIAL

    implicit none

    TREAL, dimension(isize_field), intent(IN) :: rho, u, v, w, p, e
    TREAL, dimension(isize_field, *), intent(IN) :: z1
    TREAL, dimension(isize_field), intent(OUT) :: h0, h1, h2, h3, h4
    TREAL, dimension(isize_field, *), intent(OUT) :: zh1
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5

! -------------------------------------------------------------------
    TINTEGER bcs(2, 1), i, is
    TREAL g1, g2, g3, prefactor, dummy

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

    bcs = 0

    g1 = buoyancy%vector(1)
    g2 = buoyancy%vector(2)
    g3 = buoyancy%vector(3)
    prefactor = (gama0 - C_1_R)*mach*mach

! ###################################################################
! Terms \rho u in mass, u-momentum and energy equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = C_05_R*rho(i)*u(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i) + p(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp5)
    h0 = h0 - C_2_R*tmp5

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    do i = 1, imax*jmax*kmax
        h1(i) = h1(i) - (tmp2(i) + tmp3(i) + tmp4(i) + u(i)*tmp5(i)) + g1*rho(i)
        h2(i) = h2(i) - v(i)*tmp5(i)
        h3(i) = h3(i) - w(i)*tmp5(i)
    end do

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
    if (imode_eqns == DNS_EQNS_TOTAL) then
        h4 = h4 - (e + prefactor*C_05_R*(u*u + v*v + w*w))*tmp5
! Internal energy
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        h4 = h4 - e*tmp5
    end if

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
    do is = 1, inb_scal
        zh1(:, is) = zh1(:, is) - z1(:, is)*tmp5
    end do

! ###################################################################
! Terms \rho v in mass, v-momentum and energy equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = C_05_R*rho(i)*v(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i) + p(i)
        tmp1(i) = tmp4(i)*u(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp5)
    h0 = h0 - C_2_R*tmp5

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    do i = 1, imax*jmax*kmax
        h1(i) = h1(i) - u(i)*tmp5(i)
        h2(i) = h2(i) - (tmp2(i) + tmp3(i) + tmp4(i) + v(i)*tmp5(i)) + g2*rho(i)
        h3(i) = h3(i) - w(i)*tmp5(i)
    end do

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
    if (imode_eqns == DNS_EQNS_TOTAL) then
        h4 = h4 - (e + prefactor*C_05_R*(u*u + v*v + w*w))*tmp5
! Internal energy
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        h4 = h4 - e*tmp5
    end if

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
    do is = 1, inb_scal
        zh1(:, is) = zh1(:, is) - z1(:, is)*tmp5
    end do

! ###################################################################
! Terms \rho w in mass, w-momentum and energy equations
! ###################################################################
    do i = 1, imax*jmax*kmax
        tmp4(i) = C_05_R*rho(i)*w(i)
        tmp3(i) = tmp4(i)*w(i) + p(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i)
    end do

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp5)
    h0 = h0 - C_2_R*tmp5

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    do i = 1, imax*jmax*kmax
        h1(i) = h1(i) - u(i)*tmp5(i)
        h2(i) = h2(i) - v(i)*tmp5(i)
        h3(i) = h3(i) - (tmp2(i) + tmp3(i) + tmp4(i) + w(i)*tmp5(i)) + g3*rho(i)
    end do

! -------------------------------------------------------------------
! energy equation
! -------------------------------------------------------------------
! Total energy
    if (imode_eqns == DNS_EQNS_TOTAL) then
        h4 = h4 - (e + prefactor*C_05_R*(u*u + v*v + w*w))*tmp5
! Internal energy
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        h4 = h4 - e*tmp5
    end if

! -------------------------------------------------------------------
! scalar equations
! -------------------------------------------------------------------
    do is = 1, inb_scal
        zh1(:, is) = zh1(:, is) - z1(:, is)*tmp5
    end do

! ###################################################################
! Energy equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
    if (imode_eqns == DNS_EQNS_TOTAL) then
        do i = 1, imax*jmax*kmax
            dummy = C_05_R*rho(i)*(e(i) + prefactor*C_05_R*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))) + prefactor*p(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
            tmp5(i) = (dummy - prefactor*p(i))/rho(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        h4 = h4 - (tmp2 + tmp3 + tmp4) + prefactor*rho*(g1*u + g2*v + g3*w)

! -------------------------------------------------------------------
! Internal energy formulation
! the term p div u is add in RHS_FLOW_VISCOUS to save derivatives
! -------------------------------------------------------------------
    else if (imode_eqns == DNS_EQNS_INTERNAL) then
        do i = 1, imax*jmax*kmax
            dummy = C_05_R*rho(i)*e(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
            tmp5(i) = dummy/rho(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        h4 = h4 - (tmp2 + tmp3 + tmp4)
    end if

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! ###################################################################
! energy equation (array tmp5)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp5, tmp2)
    h4 = h4 - rho*(u*tmp2 + v*tmp3 + w*tmp4)

! u-momentum equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp2)
    h1 = h1 - C_05_R*rho*(u*tmp2 + v*tmp3 + w*tmp4)

! v-momentum equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp2)
    h2 = h2 - C_05_R*rho*(u*tmp2 + v*tmp3 + w*tmp4)

! w-momentum equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp2)
    h3 = h3 - C_05_R*rho*(u*tmp2 + v*tmp3 + w*tmp4)

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_EULER_SKEWSYMMETRIC')
#endif

    return
end subroutine RHS_FLOW_EULER_SKEWSYMMETRIC
