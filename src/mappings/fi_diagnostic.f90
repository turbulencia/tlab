#include "dns_const.h"

!########################################################################
!#
!# Calculate diagnostic variables
!#
!########################################################################
subroutine FI_DIAGNOSTIC(nx, ny, nz, q, s)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use TLAB_VARS, only: imode_eqns, itransport, damkohler
    use TLAB_VARS, only: epbackground, pbackground
    use TLAB_ARRAYS, only: wrk3d
    use THERMO_VARS, only: imixture

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: q(nx*ny*nz, inb_flow_array)
    real(wp), intent(INOUT) :: s(nx*ny*nz, inb_scal_array)

    ! -------------------------------------------------------------------
    ! ###################################################################
    select case (imode_eqns)
    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then ! Calculate q_l
            call THERMO_AIRWATER_PH(nx, ny, nz, s(1, 2), s(1, 1), epbackground, pbackground)

        else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            call THERMO_AIRWATER_LINEAR(nx, ny, nz, s, s(1, inb_scal_array))

        end if

    case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
#define e(j)    q(j,4)
#define rho(j)  q(j,5)
#define p(j)    q(j,6)
#define T(j)    q(j,7)
#define vis(j)  q(j,8)

        call THERMO_CALORIC_TEMPERATURE(nx, ny, nz, s, e(1), rho(1), T(1), wrk3d)
        call THERMO_THERMAL_PRESSURE(nx, ny, nz, s, rho(1), T(1), p(1))
     if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) call THERMO_VISCOSITY(nx, ny, nz, T(1), vis(1))

    end select

    return
end subroutine FI_DIAGNOSTIC
