#include "dns_const.h"

!########################################################################
!# Skewsymmetric formulation according to Erlebacher, 1992
!# Derived from RHS_FLOW_EULER_DIVERGENCE, 3 additional derivative operations are added.
!# The mass conservation terms are implemented in the routine RHS_FLOW_EULER_SKEWSYMMETRIC
!########################################################################
subroutine RHS_SCAL_EULER_SKEWSYMMETRIC(is)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use FDM, only: g
    use TLab_Pointers
    use TLab_Arrays, only: s
    use DNS_ARRAYS, only: hs
    use OPR_PARTIAL
    implicit none

    integer, intent(in) :: is

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i
    real(wp) dummy

! ###################################################################
! divergence part
    do i = 1, imax*jmax*kmax
        dummy = 0.5_wp*rho(i)*s(i, is)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
    end do
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    hs(:, is) = hs(:, is) - (tmp2 + tmp3 + tmp4)

! convective part
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(:,is), tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:,is), tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(:,is), tmp2)
    hs(:, is) = hs(:, is) - 0.5_wp*rho*(u*tmp2 + v*tmp3 + w*tmp4)

    return
end subroutine RHS_SCAL_EULER_SKEWSYMMETRIC
