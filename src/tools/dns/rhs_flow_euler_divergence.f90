#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# 15 first derivative operations.
!#
!########################################################################
subroutine RHS_FLOW_EULER_DIVERGENCE(rho, u, v, w, p, e, h0, h1, h2, h3, h4, tmp1, tmp2, tmp3, tmp4, tmp5)

    use TLAB_VARS, only: imax, jmax, kmax, isize_field, imode_eqns
    use TLAB_VARS, only: g, buoyancy
    use TLAB_VARS, only: mach
    use THERMO_VARS, only: gama0
    use OPR_PARTIAL
    implicit none

    TREAL, dimension(isize_field), intent(IN) :: rho, u, v, w, p, e
    TREAL, dimension(isize_field), intent(OUT) :: h0, h1, h2, h3, h4
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5

! -------------------------------------------------------------------
    TINTEGER bcs(2, 1), i
    TREAL g1, g2, g3, prefactor, dummy

! ###################################################################
    bcs = 0

    g1 = buoyancy%vector(1)
    g2 = buoyancy%vector(2)
    g3 = buoyancy%vector(3)

! ###################################################################
! Terms \rho u in mass and u-momentum equations
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
    h0 = h0 - tmp5

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    h1 = h1 - (tmp2 + tmp3 + tmp4) + g1*rho

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
    h0 = h0 - tmp5

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    h2 = h2 - (tmp2 + tmp3 + tmp4) + g2*rho

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
    h0 = h0 - tmp5

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    h3 = h3 - (tmp2 + tmp3 + tmp4) + g3*rho

! ###################################################################
! Total enery equation
! ###################################################################
! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
    if (imode_eqns == DNS_EQNS_TOTAL) then
        prefactor = (gama0 - C_1_R)*mach*mach
        do i = 1, imax*jmax*kmax
            dummy = rho(i)*(e(i) + prefactor*C_05_R*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))) + prefactor*p(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
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
            dummy = rho(i)*e(i)
            tmp3(i) = dummy*w(i)
            tmp2(i) = dummy*v(i)
            tmp1(i) = dummy*u(i)
        end do
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        h4 = h4 - (tmp2 + tmp3 + tmp4)

    end if

    return
end subroutine RHS_FLOW_EULER_DIVERGENCE
