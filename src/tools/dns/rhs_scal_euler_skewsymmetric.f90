#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Skewsymmetric formulation according to Erlebacher, 1992
!# Derived from RHS_FLOW_EULER_DIVERGENCE, 3 additional derivative operations are added.
!# The mass conservation terms are implemented in the routine RHS_FLOW_EULER_SKEWSYMMETRIC
!#
!########################################################################
subroutine RHS_SCAL_EULER_SKEWSYMMETRIC(rho, u, v, w, z1, zh1, tmp1, tmp2, tmp3, tmp4)

    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use OPR_PARTIAL

    implicit none

    TREAL, dimension(isize_field), intent(IN) :: rho, u, v, w, z1
    TREAL, dimension(isize_field), intent(OUT) :: zh1
    TREAL, dimension(isize_field), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! -------------------------------------------------------------------
    TINTEGER bcs(2, 1), i
    TREAL dummy

! ###################################################################
! divergence part
    do i = 1, imax*jmax*kmax
        dummy = C_05_R*rho(i)*z1(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
    end do
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    zh1 = zh1 - (tmp2 + tmp3 + tmp4)

! convective part
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), z1, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), z1, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), z1, tmp2)
    zh1 = zh1 - C_05_R*rho*(u*tmp2 + v*tmp3 + w*tmp4)

    return
end subroutine RHS_SCAL_EULER_SKEWSYMMETRIC
