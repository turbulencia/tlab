#include "dns_const.h"

!########################################################################
!#
!# Calculate diagnostic variables
!#
!########################################################################
subroutine FI_DIAGNOSTIC(nx, ny, nz, q, s)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use TLAB_VARS, only: imode_eqns, itransport, damkohler, buoyancy
    use TLAB_ARRAYS, only: wrk1d, wrk3d
    use FI_SOURCES
    use Thermodynamics, only: imixture
    use THERMO_THERMAL
    use THERMO_CALORIC
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use AVGS, only: AVG1V2D_V

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: q(nx*ny*nz, inb_flow_array)
    real(wp), intent(INOUT) :: s(nx*ny*nz, inb_scal_array)

    ! -------------------------------------------------------------------
    ! ###################################################################
    select case (imode_eqns)
    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        ! Calculate liquid content q_l
        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
            call THERMO_ANELASTIC_PH(nx, ny, nz, s(:, 2), s(1, 1))

        else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            call THERMO_AIRWATER_LINEAR(nx*ny*nz, s, s(:, inb_scal_array))

        end if

        ! Calculate mean bbackground
        if (buoyancy%type == EQNS_BOD_NORMALIZEDMEAN .or. buoyancy%type == EQNS_BOD_SUBTRACTMEAN) then
            call AVG1V2D_V(nx, ny, nz, 1, s(:, 1), bbackground, wrk1d)
        end if

    case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
#define e(j)    q(j,4)
#define rho(j)  q(j,5)
#define p(j)    q(j,6)
#define T(j)    q(j,7)
#define vis(j)  q(j,8)

        call THERMO_CALORIC_TEMPERATURE(nx*ny*nz, s, e(1), rho(1), T(1), wrk3d)
        call THERMO_THERMAL_PRESSURE(nx*ny*nz, s, rho(1), T(1), p(1))
        if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) call THERMO_VISCOSITY(nx*ny*nz, T(1), vis(1))

    end select

    return
end subroutine FI_DIAGNOSTIC
